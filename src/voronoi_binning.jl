
abstract type BinningStrategy end
abstract type VoronoiStrategy <: BinningStrategy end

struct InitialVoronoi <: VoronoiStrategy end
struct CentroidalVoronoi <: VoronoiStrategy end
struct WeightedVoronoi <: VoronoiStrategy 
    weight_function::Base.Callable
    function WeightedVoronoi(weight_function=default_WVT_func)
        new(weight_function)
    end
end

voronoi2Dbinning(
    x::AbstractVector{T1},
    y::AbstractVector{T1},
    S::AbstractVector{T2},
    N::AbstractVector{T2},
    target_SN::T2;
    kwargs...) where {
        T1<:Real,
        T2<:Real
    } = voronoi2Dbinning(x, y, S, N, target_SN, default_SN_func, WeightedVoronoi(); kwargs...)

voronoi2Dbinning(
    x::AbstractVector{T1},
    y::AbstractVector{T1},
    S::AbstractVector{T2},
    N::AbstractVector{T2},
    target_SN::T2,
    bin_strategy::VoronoiStrategy;
    kwargs...) where {
        T1<:Real,
        T2<:Real
    } = voronoi2Dbinning(x, y, S, N, target_SN, default_SN_func, bin_strategy; kwargs...)

voronoi2Dbinning(
    x::AbstractVector{T1},
    y::AbstractVector{T1},
    S::AbstractVector{T2},
    N::AbstractVector{T2},
    target_SN::T2,
    SN_func::Base.Callable;
    kwargs...) where {
        T1<:Real,
        T2<:Real
    } = voronoi2Dbinning(x, y, S, N, target_SN, SN_func, WeightedVoronoi(); kwargs...)

voronoi2Dbinning(
    x::AbstractVector{T1}, 
    y::AbstractVector{T1}, 
    S::AbstractVector{T2}, 
    N::AbstractVector{T2}, 
    target_SN::T2,
    pixel_size::T3;
    kwargs...) where {
        T1<:Real, 
        T2<:Real,
        T3<:Real
    } = voronoi2Dbinning(x, y, S, N, target_SN, default_SN_func, pixel_size, WeightedVoronoi(); kwargs...)

voronoi2Dbinning(
    x::AbstractVector{T1}, 
    y::AbstractVector{T1}, 
    S::AbstractVector{T2}, 
    N::AbstractVector{T2}, 
    target_SN::T2,
    pixel_size::T3,
    bin_strategy::VoronoiStrategy;
    kwargs...) where {
        T1<:Real, 
        T2<:Real,
        T3<:Real
    } = voronoi2Dbinning(x, y, S, N, target_SN, default_SN_func, pixel_size, bin_strategy; kwargs...)

voronoi2Dbinning(
    x::AbstractVector{T1}, 
    y::AbstractVector{T1}, 
    S::AbstractVector{T2}, 
    N::AbstractVector{T2}, 
    target_SN::T2,
    SN_func::Base.Callable,
    pixel_size::T3;
    kwargs...) where {
        T1<:Real, 
        T2<:Real,
        T3<:Real
    } = voronoi2Dbinning(x, y, S, N, target_SN, SN_func, pixel_size, WeightedVoronoi(); kwargs...)

function voronoi2Dbinning(
    x::AbstractVector{T1}, 
    y::AbstractVector{T1}, 
    S::AbstractVector{T2}, 
    N::AbstractVector{T2}, 
    target_SN::T2,
    SN_func::Base.Callable,
    bin_strategy::VoronoiStrategy;
    kwargs...) where {
        T1<:Real, 
        T2<:Real
    }
    if length(x) > 10^4
        error("Dataset has more than 10^4 pixels. Please provide a \'pixel_size\' argument.")
    end

    pixel_size = Inf
    for i in 1:length(x)
        for j in (i+1):length(x)
            sep = hypot(x[j]-x[i], y[j]-y[i])
            pixel_size = sep < pixel_size ? sep : pixel_size
        end
    end

    voronoi2Dbinning(x, y, S, N, target_SN, SN_func, pixel_size, bin_strategy; kwargs...)
end

function voronoi2Dbinning(
    x::AbstractVector{T1}, 
    y::AbstractVector{T1}, 
    S::AbstractVector{T2}, 
    N::AbstractVector{T2}, 
    target_SN::T2,
    SN_func::Base.Callable,
    pixel_size::T3,
    bin_strategy::InitialVoronoi;
    verbose::Bool=false,
    plot::Bool=false) where {
        T1<:Real, 
        T2<:Real,
        T3<:Real
    }
    t1 = time()
    bin_numbers, xnode, ynode = initial_voronoi_step(x, y, S, N, target_SN, SN_func, pixel_size; verbose=verbose)
    t2 = time()
    weights = ones(eltype(xnode), length(xnode))
    bin_numbers, x̄, ȳ, SN, area = get_bin_quantities(x, y, S, N, xnode, ynode, weights, SN_func, pixel_size)
    single = area .≈ pixel_size.^2
    if verbose
        println("Unbinned pixels: $(sum(single)) / $(length(x))")
        println("Fractional S/N scatter: ", std(SN[.~single] .- target_SN)/target_SN * 100, " %")
        println("Elapsed time for accretion: $(@sprintf "%.2f" (t2 - t1)) seconds")
    end
    if plot
        if !DO_PLOTTING
            println("WARNING: plot=true argument was passed, but the PyPlot package was not found - plotting will be skipped!")
            fig, ax = nothing, nothing
        else
            fig, ax = plot_voronoi_tessellation(x, y, S, N, bin_numbers, xnode, ynode, x̄, ȳ, area, SN, target_SN, pixel_size)
        end
        return bin_numbers, xnode, ynode, x̄, ȳ, SN, area, weights, fig, ax
    end
    bin_numbers, xnode, ynode, x̄, ȳ, SN, area, weights
end


"""
    voronoi2Dbinning(x, y, S, N, target_SN[, SN_func, pixel_size, bin_strategy]; verbose=false, plot=false)

Calculate a Voronoi binning scheme for 2-dimensional data such that the signal to noise ratio of each bin is roughly constant.
This function implements the adaptive Voronoi binning schemes described in Cappellari & Copin (2003) and Diehl & Statler (2006),
while allowing configuration options to modify how the data is binned.

# Arguments

- `x::AbstractVector{T1}`: The x-coordinates of each pixel to be binned, flattened to a 1-dimensional vector. The pixel grid
    is assumed to be regular, but can contain holes or an irregular boundary. The units are arbitrary and can be anything, as long
    as they are consistent with the units of `y`.
- `y::AbstractVector{T1}`: The y-coordinates of each pixel to be binned, flattened to a 1-dimensional vector. Like `x`, 
    the units are arbitrary. Must be the same length as `x`.
- `S::AbstractVector{T2}`: The signal of each pixel that has coordinates (x, y). Must be the same length as `x` and `y`, such that
    each element of `S` corresponds to an element of `x` and `y`. This is used to calculate the signal-to-noise ratios.
- `N::AbstractVector{T2}`: The noise of each pixel that has coordinates (x, y). Must be the same length as `x` and `y`, such that
    each element of `N` corresponds to an element of `x` and `y`. This is used to calculate the signal-to-noise ratios.
- `target_SN::T2`: The target signal-to-noise ratio that one wishes to achieve during the binning process. Individual pixels that have
    a higher S/N than this will not be binned, while pixels with S/Ns lower than this will be binned to their neighbors until the bin
    has a sufficiently close S/N.

# Optional Arguments

- `SN_func::Base.Callable`: A callable function that computes the signal-to-noise ratio of a set of pixels. The arguments to the function
    must be `(index, S, N)` where `index` is a set of indices for `S` and `N` that describes a bin or partial bin. If not specified, the
    default function calculates the S/N of a set of pixels with signals (S1, S2, ...) and noises (N1, N2, ...) as 
    (S1 + S2 + ...) / √(N1^2 + N2^2 + ...). However, the quantity returned by the SN_func need not necessarily be a S/N ratio. It can,
    in principle, be any quantity that the user wishes to equalize over the voronoi bins.
- `pixel_size::T3`: A real number giving the size of each pixel. If the units are pixels, this should be 1.0. If not provided, it will
    be calculated automatically from the `x` and `y` vectors, assuming they are small enough for the computation time to be reasonable.
- `bin_strategy::VoronoiStrategy`: A struct indicating the type of binning strategy to use. Can be InitialVoronoi, CentroidalVoronoi, or
    WeightedVoronoi. InitialVoronoi simply computes initial voronoi bins using the accretion algorithm specified by Cappellari & Copin (2003).
    CentroidalVoronoi takes these initial bins and performs an iterative adjustment procedure using a modified Lloyd algorithm with the
    S/N density squared. WeightedVoronoi similarly uses an iterative adjustment procedure on the initial bins, but uses a weighting scheme,
    which by default is the weighting scheme proposed by Diehl & Statler (2006), but can be modified to be any arbitrary weighting scheme by
    specifying the weight function with WeightedVoronoi(weight_function). The weight function must take arguments (indices, S, N, SN_func). 
    The default is WeightedVoronoi with the weighting scheme of Diehl & Statler (2006).

# Keyword Arguments

- `verbose::Bool=false`: A flag that enables progress and performance messages.
- `dist2_thresh::T3=1.44pixel_size^2`: Threshold on the distance^2 to consider during the bin accretion step. If a new pixel's distance^2 is within
    this threshold, then it may be accreted.
- `R_thresh::T3=0.3`: Threshold on the roundness parameter to consider during the bin accretion step. If a new pixel does not increase the roundness
    above R_thresh, then it may be accreted.
- `cvt_tolerance::T2=0.1`: The tolerance level to check against during the iterative bin adjustment procedure, if the `bin_strategy` is CentroidalVoronoi
    or WeightedVoronoi. This is a tolerance on the sum of the differences between the positions of the node generators at the current step and the previous step.
- `plot::Bool=false`: A flag that, if set, will automatically plot the voronoi bin map as well as the S/N of each bin as a function of radius.

`{T1<:Integer, T2<:Real, T3<:Real}`

# Returns

- `bin_numbers::Vector{Int}`: This vector is the same length as `x` and `y` and contains integers that act as unique identifiers for each bin.
    If a pixel located at (x, y) has index `i` in the `x` and `y` arrays, then it is located in the bin labeled by `bin_numbers[i]`. This is
    all one needs to fully specify the binned data, but other outputs are provided for utility and completeness.
- `xnode::Vector{<:Real}`: The vector of x positions for the generator points for each bin, which can be used with the voronoi_tessellation
    function to generate the bins_numbers. N.B. these are *not* the same as the bin centroids.
- `ynode::Vector{<:Real}`: Like `xnode`, but for the y positions.
- `x̄::Vector{<:Real}`: The vector of x positions for the centroids (weighted means) of the final set of bins.
- `ȳ::Vector{<:Real}`: Like `x̄`, but for the y positions.
- `SN::Vector{<:Real}`: The S/N values (or whatever values are returned by the SN_func) for the final set of bins.
- `area::Vector{<:Real}`: The areas of each bin, in units given by the units of `x` and `y` squared.
- `weights::Vector{<:Real}`: If using the WeightedVoronoi strategy, this gives the weights used to compute the `xnode` and `ynode` positions.
- `[fig::PyObject]`: If `plot=true`, the figure and axes matplotlib objects of the plot are returned so they can be further modified.
- `[ax::PyObject]`: If `plot=true`, the figure and axes matplotlib objects of the plot are returned so they can be further modified.


# Notes from the original Python vorbin package:

    1. Be careful when attempting to bin pixels with almost no signal but considerable noise. In these situations, it is recommended
       to make a S/N cut and only bin the pixels that satisfy this cut -- however, if this creates pixels in disconnected regions, each
       region must be binned separately. Otherwise, one may try to optimally weight the pixels before binning -- refer to section 2.1 of
       Cappellari & Copin (2003).

    2. Ignore (#1) for X-ray or other data where photons are counted individually and Poissonian statistics is appropriate. This is
       because the S/N of a bin can never decrease by adding a pixel, so it is preferrable to include the pixels with 0 signal.

    3. For very large images where computation time becomes an issue, it may be preferrable to approach the binning hierarchically. 
       This would be done by first rebinning the image regularly (i.e. by a factor of 4x4) and applying the voronoi binning scheme on
       the binned image. Then, take all of the unbinned pixels of the voronoi-binned result and transform them back into their original 
       individual, full-resolution pixels, and apply the voronoi binning on each connected region of pixels.

"""
function voronoi2Dbinning(
    x::AbstractVector{T1}, 
    y::AbstractVector{T1}, 
    S::AbstractVector{T2}, 
    N::AbstractVector{T2}, 
    target_SN::T2,
    SN_func::Base.Callable,
    pixel_size::T3,
    bin_strategy::Union{CentroidalVoronoi,WeightedVoronoi};
    verbose::Bool=false,
    dist2_thresh::T3=(1.2pixel_size)^2,
    R_thresh::T3=0.3,
    cvt_tolerance::T2=0.1,
    plot::Bool=false) where {
        T1<:Real, 
        T2<:Real,
        T3<:Real
    }
    t1 = time()
    bin_numbers, xnode, ynode = initial_voronoi_step(x, y, S, N, target_SN, SN_func, pixel_size; verbose=verbose, dist2_thresh=dist2_thresh,
        R_thresh=R_thresh)
    t2 = time()
    verbose && println("Performing modified Lloyd algorithm:")
    xnode, ynode, weights, iters = voronoi_bin_adjustment!(x, y, S, N, xnode, ynode, pixel_size, SN_func, bin_strategy, cvt_tolerance; verbose=verbose)
    t3 = time()
    verbose && println("Finished in $iters iterations.")
    bin_numbers, x̄, ȳ, SN, area = get_bin_quantities(x, y, S, N, xnode, ynode, weights, SN_func, pixel_size)
    single = area .≈ pixel_size.^2
    if verbose
        println("Unbinned pixels: $(sum(single)) / $(length(x))")
        println("Fractional S/N scatter: ", std(SN[.~single] .- target_SN)/target_SN * 100, " %")
        println("Elapsed time for accretion: $(@sprintf "%.2f" (t2 - t1)) seconds")
        println("Elapsed time for optimization: $(@sprintf "%.2f" (t3 - t2)) seconds")
    end
    if plot
        if !DO_PLOTTING
            fig, ax = nothing, nothing
        else
            fig, ax = plot_voronoi_tessellation(x, y, S, N, bin_numbers, xnode, ynode, x̄, ȳ, area, SN, target_SN, pixel_size)
        end
        return bin_numbers, xnode, ynode, x̄, ȳ, SN, area, weights, fig, ax
    end
    bin_numbers, xnode, ynode, x̄, ȳ, SN, area, weights
end


