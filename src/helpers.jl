
"""
    check_inputs(x, y, S, N, target_SN, SN_func)

Helper function to perform a few checks to make sure the inputs make sense. Namely, checking that
x, y, S, and N are all the same length, and making sure there is enough S/N to meet the target_SN
requirement.
"""
function check_inputs(
    x::AbstractVector{T1},
    y::AbstractVector{T1},
    S::AbstractVector{T2},
    N::AbstractVector{T2}, 
    target_SN::T2,
    SN_func::Base.Callable) where {
        T1<:Integer,
        T2<:Real
    }

    # Input checking
    @assert length(x) == length(y) == length(S) == length(N) "x, y, S, and N must all have the same length!"
    @assert all((N .> 0) .& isfinite.(N)) "Noise values must all be finite and positive!"

    # Making sure we can reach the target S/N
    SN_func(findall(N .> 0), S, N) < target_SN && error("There is not enough S/N in the entire set of pixels to satisfy" * 
        " the provided target S/N. It is possible that many pixels have noise but little to no signal. These pixels should" *
        " either be removed or optimally weighted -- refer to Cappellari & Copin (2003, Sec.2.1).")
    # Making sure that binning is actually necessary
    minimum(S./N) > target_SN && error("All pixels meet the target S/N threshold. No binning is required.")

end


"""
    default_SN_func(indices, S, N)

The default function for calculating the S/N of multiple pixels, using the standard definition
of (S1 + S2 + ...) / √(N1^2 + N2^2 + ...)
"""
@inline @views default_SN_func(
    indices::AbstractVector{T1}, 
    S::AbstractVector{T2}, 
    N::AbstractVector{T2}) where {
        T1<:Integer,
        T2<:Real
    } = sum(S[indices]) / √(sum(x -> x^2, N[indices]))


"""
    default_WVT_func(indices, S, N, SN_func)

The default function for calculating the weights during the weighted voronoi tessellation step.
The weights are those proposed by Diehl & Statler (2006), i.e. (S/N) / area

N.B. `indices` is expected to be a BITVECTOR giving the positions within `S` and `N`.
"""
@inline default_WVT_func(
    indices::BitVector,
    S::AbstractVector{T2},
    N::AbstractVector{T2},
    SN_func::Base.Callable) where {
        T2<:Real
    } = SN_func(indices, S, N) / sum(indices)


"""
    voronoi_tessellation(x, y, xnode, ynode[, weights])

Computes the (optionally weighted) Voronoi tessellation of a grid of pixels
given node points.
"""
function voronoi_tessellation(
    x::AbstractVector{T1},
    y::AbstractVector{T1},
    xnode::AbstractVector{T2},
    ynode::AbstractVector{T2}) where {
        T1<:Real,
        T2<:Real
    }
    # N.B. minor differences between NearestNeighbor.jl's KDTree and scipy's cKDTree cause the 
    # binning results to be slightly different from the python vorbin package, but the S/N scatter
    # remains consistent within the percent level.
    tree = KDTree([xnode ynode]'; leafsize=16)  # use leafsize=16 to match the default for scipy's KDTree
    points = [x y]'
    nn(tree, points)[1]

    # using cKDTree from scipy reproduces python vorbin exactly:
    # tree = py_spatial.cKDTree([xnode ynode])
    # tree.query([x y])[2] .+ 1
end

function voronoi_tessellation(
    x::AbstractVector{T1},
    y::AbstractVector{T1},
    xnode::AbstractVector{T2},
    ynode::AbstractVector{T2},
    weights::AbstractVector{T3}) where {
        T1<:Real,
        T2<:Real,
        T3<:Real
    } 
    if all(isone.(weights))
        return voronoi_tessellation(x, y, xnode, ynode)
    end
    indices = zeros(Int64, length(x))
    @inbounds @simd for j in eachindex(x)
        xj, yj = x[j], y[j]
        minind = 1
        minval = ((xj - xnode[1])^2 + (yj - ynode[1])^2) * weights[1]
        for k ∈ eachindex(xnode)
            dist = ((xj - xnode[k])^2 + (yj - ynode[k])^2) * weights[k]
            newmin = dist < minval
            minval = newmin ? dist : minval
            minind = newmin ? k : minind
        end
        indices[j] = minind
        # indices[j] = argmin(@. ((x[j] - xnode)^2 + (y[j] - ynode)^2) * weights)
    end
    indices
end


"""
    roundness(x, y, pixel_size)

Compute the roundness metric, R = R_max/R_eff - 1, of a bin, which is a measure of how
close the bin shape is to a circle, in which case R_max = R_eff and R = 0.

See equation (5) of Cappellari & Copin (2003)
"""
function roundness(
    x::AbstractVector{T1},
    y::AbstractVector{T1},
    pixel_size::T2) where {
        T1<:Integer,
        T2<:Real
    }
    r_eff = √(length(x)/π)*pixel_size   # effective radius of a disc with the same area as the bin
    x̄, ȳ = mean(x), mean(y)

    # r_max = √(maximum((x .- x̄).^2 .+ (y .- ȳ).^2))
    # make this fast and more memory efficient with LoopVectorization
    r_max = nextfloat(typemin(Float64))
    @turbo for i in eachindex(x)
        r_test = (x[i] - x̄)^2 + (y[i] - ȳ)^2
        r_max = r_test > r_max ? r_test : r_max
    end
    r_max = √(r_max)

    r_max/r_eff - 1
end


"""
    mindist_maximumoverdrive(x, y, x₀, y₀)

A fast re-implementation of `findmin(@. (x - x₀)^2 + (y - y₀)^2)` with almost no allocations
"""
function mindist_maximumoverdrive(x, y, x₀, y₀)
    mindist = (x[1] - x₀)^2 + (y[1] - y₀)^2
    minind = 1
    @inbounds for i ∈ eachindex(x)
        dist = (x[i] - x₀)^2 + (y[i] - y₀)^2
        newmin = dist < mindist
        mindist = newmin ? dist : mindist
        minind = newmin ? i : minind
    end
    return mindist, minind
end


"""
    voronoi_bin_accretion(x, y, S, N, target_SN, SN_func, pixel_size; verbose=false)

Computes the initial voronoi bins using an accretion process described in section 5.1 of 
Cappellari & Copin (2003).  If the binning strategy is InitialVoronoi(), these are the final 
output Voronoi bins.
"""
function voronoi_bin_accretion(
    x::AbstractVector{T1},
    y::AbstractVector{T1},
    S::AbstractVector{T2},
    N::AbstractVector{T2}, 
    target_SN::T2,
    SN_func::Base.Callable,
    pixel_size::T3;
    verbose::Bool=false,
    dist2_thresh::Real=1.44pixel_size^2,
    R_thresh::Real=0.3) where {
        T1<:Integer,
        T2<:Real,
        T3<:Real
    } 
    # Prepare outputs
    n = length(x)
    bin_numbers = zeros(Int, n)            # bin number associated with each (x,y) point
    flags = falses(n)                        # flag for bins that have been marked as good 

    # Rough estimate of the expected final number of bins
    if verbose
        w = (S ./ N) .< target_SN
        max_N_guess = round(Int, sum((S[w]./N[w]).^2) / target_SN^2 + sum(.~w))
    end

    # (i) start the first bin from the highest S/N pixel of the image
    current_bin = Int64[argmax(S ./ N)]
    current_SN = SN_func(current_bin, S, N)

    # Prepare arrays for the currently unbinned and binned pixels
    unbinned = findall(iszero, bin_numbers) 
    binned = [current_bin[1]]
    deleteat!(unbinned, current_bin[1])

    function test_add_pixel(bin::Vector{Int}, pixel::Int, old_SN::Real)
        test_bin = [bin; pixel]
        # min_dist = @views minimum((x[bin] .- x[pixel]).^2 .+ (y[bin] .- y[pixel]).^2)
        # LoopVectorization to the rescue
        min_dist = prevfloat(typemax(Float64))
        @turbo for bi ∈ eachindex(bin)
            dist_test = (x[bin[bi]] - x[pixel])^2 + (y[bin[bi]] - y[pixel])^2
            min_dist = dist_test < min_dist ? dist_test : min_dist
        end
            
        R = roundness(x[test_bin], y[test_bin], pixel_size)
        new_SN = SN_func(test_bin, S, N)

        # Continue the accretion process if the following criteria all hold:
        #     (a) TOPOLOGICAL: the new pixel is adjacent to the current bin
        #     (b) MORPHOLOGICAL: adding the new pixel keeps the roundness below a given threshold
        #     (c) UNIFORMITY: adding the new pixel pushes the S/N closer to the target S/N 
        topology = min_dist > dist2_thresh
        morphology = R > R_thresh
        uniformity = (abs(new_SN - target_SN) > abs(old_SN - target_SN)) || (old_SN > new_SN)
        
        # Return true if the pixel should be accreted, also return the new S/N should the pixel be accreted
        !(topology || morphology || uniformity), new_SN
    end

    ispositive = x -> x > 0

    # Assign the first bin to bin_number = 1
    # With N pixels, there will be no more than N bins
    for bin_index ∈ 1:n

        verbose && print("  $bin_index / $max_N_guess \r")

        bin_numbers[current_bin] .= bin_index   # here current_bin only contains 1 pixel
        x̄ = x[current_bin[1]]
        ȳ = y[current_bin[1]]  # the "centroids" are just the pixel coordinates

        while true
            # Stop the loop if we have binned all pixels
            if all(ispositive, bin_numbers)
                break
            end

            # (ii) Find the closest unbinned pixel 
            # candidate = @views argmin((x[unbinned] .- x̄).^2 .+ (y[unbinned] .- ȳ).^2)
            x_unbinned = @view x[unbinned]
            y_unbinned = @view y[unbinned]
            _, candidate = mindist_maximumoverdrive(x_unbinned, y_unbinned, x̄, ȳ)
            test_pixel = unbinned[candidate]

            # (iii) Test for acceptance
            last_SN = current_SN
            add_pixel, current_SN = test_add_pixel(current_bin, test_pixel, last_SN)

            if !add_pixel
                # (iv) the accretion process of the current bin ends if the bin S/N is greater than a given threshold (e.g. 80%)
                # of the target S/N
                if last_SN > 0.8*target_SN
                    flags[current_bin] .= 1
                end
                break
            end

            # Otherwise we want to add the candidate pixel to the current bin and continue the accretion
            bin_numbers[test_pixel] = bin_index
            current_bin = [current_bin; test_pixel]
            # Update the unbinned and binned arrays
            deleteat!(unbinned, candidate)
            push!(binned, test_pixel)

            # Update the centroid 
            @views x̄, ȳ = mean(x[current_bin]), mean(y[current_bin])
        end

        # Stop if all pixels have been binned
        if length(binned) == n
            break
        end

        # (v) evaluate the mass centroid of all the binned pixels
        @views x̄, ȳ = mean(x[binned]), mean(y[binned])

        # Find the closest unbinned pixel to the centroid and start a new bin from there
        # unless the remaining pixels do not have enough S/N, then stop
        if SN_func(unbinned, S, N) < target_SN
            break
        end
        # min_dist = @views argmin((x[unbinned] .- x̄).^2 .+ (y[unbinned] .- ȳ).^2)
        x_unbinned = @view x[unbinned]
        y_unbinned = @view y[unbinned]
        _, min_dist = mindist_maximumoverdrive(x_unbinned, y_unbinned, x̄, ȳ)
        current_bin = [unbinned[min_dist]]
        current_SN = SN_func(current_bin, S, N)
        deleteat!(unbinned, min_dist)
        push!(binned, current_bin[1])

    end

    verbose && println()

    # Set all bins that did not reach the target S/N to 0
    bin_numbers[.~flags] .= 0

    bin_numbers 

end


"""
    group(func, x[, groups, index])

Apply a function `func` onto a array `x` after subdividing it into groups, where `groups` gives a unique integer
label to each group and is the same shape as `x`. Optionally, `index` specifies which groups to return (giving the labels
of the groups that should be returned), otherwise all groups are returned.
"""
group(func::Base.Callable, x) = func(x)
group(func::Base.Callable, x::AbstractArray, groups::AbstractArray{<:Integer}) = group(func, x, groups, sort(unique(groups)))
function group(
    func::Base.Callable, 
    x::AbstractArray, 
    groups::AbstractArray{<:Integer}, 
    index::AbstractVector{<:Integer})

    out = zeros(length(index))
    @inbounds for (i, gi) in enumerate(index)
        w = groups .== gi
        out[i] = func(x[w])
    end
    out

end


"""
    voronoi_reassign_bad_bins!(bin_numbers, x, y)

Checks for any bin_numbers that are 0 and reassigns them to the closest bin. This implements
steps (vi) and (vii) of the initial accretion process described in section 5.1 of Cappellari and
Copin (2003).

NOTE: bin_numbers is modified in-place
"""
function voronoi_reassign_bad_bins!(
    bin_numbers::AbstractVector{T1},
    x::AbstractVector{T2},
    y::AbstractVector{T2}) where {
        T1<:Integer,
        T2<:Integer
    }
    # Get the centroids of all the successful bins
    good = sort(unique(bin_numbers[bin_numbers .> 0]))
    xnode = group(mean, x, bin_numbers, good)
    ynode = group(mean, y, bin_numbers, good)

    # (vi) reassign the unsuccessfully binned pixels to the closest centroid
    bad = iszero.(bin_numbers)
    indices = voronoi_tessellation(x[bad], y[bad], xnode, ynode)
    bin_numbers[bad] = good[indices]

    # (vii) recompute the centroids of each bin that are now different due to the previous step
    good = sort(unique(bin_numbers))
    xnode = group(mean, x, bin_numbers, good)
    ynode = group(mean, y, bin_numbers, good)

    xnode, ynode
end


"""
    voronoi_bin_adjustment!(x, y, S, N, xnode, ynode, pixel_size, SN_func, bin_strategy; verbose=false)

Further optimizes the Voronoi bins by iteratively moving the nodes. The way the nodes are moved depends on
the bin_strategy:
    - CentroidalVoronoi(): Equalize the quantity ρ² = (S/N)⁴ using a modified Lloyd algorithm.
                           See section 4.1 of Cappellari & Copin (2003).
    - WeightedVoronoi():   Equalize using the weights (S/N)/A where A is the bin area. 
                           See Diehl & Statler (2006).
"""
function voronoi_bin_adjustment!(
    x::AbstractVector{T1},
    y::AbstractVector{T1},
    S::AbstractVector{T2},
    N::AbstractVector{T2},
    xnode::AbstractVector{T3},
    ynode::AbstractVector{T3},
    pixel_size::T4,
    SN_func::Base.Callable,
    bin_strategy::VoronoiStrategy,
    tolerance::T2;
    verbose::Bool=false) where {
        T1<:Integer,
        T2<:Real,
        T3<:Real,
        T4<:Real
    }
    # Prepare the density distribution
    ρ² = (S ./ N).^4
    # Start with all scales = 1
    weights = ones(eltype(xnode), size(xnode))
    # Store differences to check for convergence
    Δ = eltype(xnode)[]
    # tol = 1e-4 * max(x |> extrema |> collect |> diff |> first |> abs, y |> extrema |> collect |> diff |> first |> abs)
    good = nothing

    for i ∈ 1:length(xnode)
        # (i) start with our initial set of generators from the bin_accretion_step
        xnode_prev, ynode_prev = copy(xnode), copy(ynode)

        # (ii) perform a voronoi tessellation with the generators
        # this will contain *on average* constant-mass bins, but with large scatter
        bin_numbers = voronoi_tessellation(x, y, xnode, ynode, weights)

        # (iii) compute the mass centroids of each bin using the density squared
        # xcentroid = ∫xρ²dr / ∫ρ²dr, ycentroid = ∫yρ²dr / ∫ρ²dr
        good = bin_numbers |> unique |> sort
        if bin_strategy isa WeightedVoronoi
            # Use modified weights as proposed by Diehl & Statler (2006)
            @inbounds for gi in good
                index = bin_numbers .== gi
                xnode[gi] = mean(x[index])
                ynode[gi] = mean(y[index])
                weights[gi] = bin_strategy.weight_function(index, S, N, SN_func)
            end
        elseif bin_strategy isa CentroidalVoronoi
            mass = group(sum, ρ², bin_numbers, good)
            xnode = group(sum, x .* ρ², bin_numbers, good) ./ mass
            ynode = group(sum, y .* ρ², bin_numbers, good) ./ mass
        else
            error("Unrecognized bin strategy: $bin_strategy")
        end

        # (iv) check if the positions of the generators have converged within a reasonable margin
        # Δi = sqrt(sum((xnode .- xnode_prev).^2 .+ (ynode .- ynode_prev).^2)) / pixel_size
        # Fast LoopVectorization implementation
        Δi = 0.
        @turbo for j ∈ eachindex(xnode)
            Δi += (xnode[j] - xnode_prev[j])^2 + (ynode[j] - ynode_prev[j])^2
        end
        Δi = √(Δi/pixel_size)
        push!(Δ, Δi)
        verbose && print("Iteration $i Difference: $(@sprintf "%.4g" Δ[i])     \r")

        # Also check if Δ repeats to avoid infinite cycling
        if Δ[i] < tolerance || Δ[i] ∈ Δ[1:i-1]
            break
        end
    end

    verbose && println()

    # If the final difference is still greater than 0, re-compute the voronoi tessellation
    if Δ[end] > 0
        bin_numbers = voronoi_tessellation(x, y, xnode, ynode, weights) 
        good = sort(unique(bin_numbers))
    end

    xnode[good], ynode[good], weights[good], length(Δ)
end


"""
    initial_voronoi_step(x, y, S, N, target_SN, SN_func, pixel_size; verbose=false)

Performs the initial bin accretion process and reassigning bad bins.
"""
function initial_voronoi_step(
    x::AbstractVector{T1}, 
    y::AbstractVector{T1}, 
    S::AbstractVector{T2}, 
    N::AbstractVector{T2}, 
    target_SN::T2,
    SN_func::Base.Callable,
    pixel_size::T3;
    kwargs...) where {
        T1<:Integer, 
        T2<:Real,
        T3<:Real
    }
    check_inputs(x, y, S, N, target_SN, SN_func)
    verbose = haskey(kwargs, :verbose) ? kwargs[:verbose] : false

    verbose && println("Bin accretion:")
    bin_numbers = voronoi_bin_accretion(x, y, S, N, target_SN, SN_func, pixel_size; kwargs...)
    verbose && println("Obtained $(maximum(bin_numbers)) initial bins.")

    verbose && println("Reassigning bad bins:")
    xnode, ynode = voronoi_reassign_bad_bins!(bin_numbers, x, y)
    verbose && println("Obtained $(length(xnode)) good bins.")

    bin_numbers, xnode, ynode
end


"""
    get_bin_quantities(x, y, S, N, xnode, ynode, weights, SN_func)

Calculate the final bin numbers, centroids, S/N, and areas.
"""
function get_bin_quantities(
    x::AbstractVector{T1}, 
    y::AbstractVector{T1}, 
    S::AbstractVector{T2}, 
    N::AbstractVector{T2}, 
    xnode::AbstractVector{T3},
    ynode::AbstractVector{T3},
    weights::AbstractVector{T3},
    SN_func::Base.Callable,
    pixel_size::T4) where {
        T1<:Integer, 
        T2<:Real,
        T3<:Real,
        T4<:Real
    }
    # Get the final bin numbers for each pixel
    bin_numbers = voronoi_tessellation(x, y, xnode, ynode, weights)
    # Compute the flux-weighted centroids and S/N of each bin
    good = sort(unique(bin_numbers))
    x̄ = group(mean, x, bin_numbers, good)
    ȳ = group(mean, y, bin_numbers, good)
    area = counts(bin_numbers) .* pixel_size.^2
    SN = zeros(eltype(xnode), length(xnode))
    good = sort(unique(bin_numbers))
    for gi in good
        w = findall(bin_numbers .== gi)
        SN[gi] = SN_func(w, S, N)
    end

    bin_numbers, x̄, ȳ, SN, area
end
