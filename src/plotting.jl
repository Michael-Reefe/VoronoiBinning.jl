

"""
    plot_voronoi_tessellation(x, y, S, N, bin_numbers, xnode, ynode, x̄, ȳ, area, SN, target_SN, pixel_size)

Make a plot of the voronoi bin map as well as the S/N of each pixel as a function of radial distance from 
the brightest bin.
"""
function plot_voronoi_tessellation(
    x::AbstractVector{T1},
    y::AbstractVector{T1},
    S::AbstractVector{T2},
    N::AbstractVector{T2},
    bin_numbers::AbstractVector{T3},
    xnode::AbstractVector{T4},
    ynode::AbstractVector{T4},
    x̄::AbstractVector{T5},
    ȳ::AbstractVector{T5},
    area::AbstractVector{<:Real},
    SN::AbstractVector{<:Real},
    target_SN::Real,
    pixel_size::Real) where {
        T1<:Integer,
        T2<:Real,
        T3<:Integer,
        T4<:Real,
        T5<:Real
    }
    labels = sortperm(rand(length(x)))[bin_numbers]
    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    nx = round(Int, (xmax - xmin)/pixel_size + 1)
    ny = round(Int, (ymax - ymin)/pixel_size + 1)
    img = ones(nx, ny) .* NaN 
    i = round.(Int, (x .- xmin)./pixel_size .+ 1)
    j = round.(Int, (y .- ymin)./pixel_size .+ 1)
    for (ii, jj, ll) in zip(i, j, labels)
        img[ii, jj] = ll
    end

    fig, ax = plt.subplots(1, 2, figsize=(12,6))
    cmap = plt.matplotlib.cm.twilight_shifted
    cmap.set_bad(color="k")
    ax[1].imshow(img', origin="lower", interpolation="nearest", cmap=cmap, 
        extent=[xmin - pixel_size/2, xmax + pixel_size/2, ymin - pixel_size/2, ymax + pixel_size/2])
    ax[1].plot(xnode, ynode, "w+", scalex=false, scaley=false)
    # ax[1].axis(:off)
    ax[1].set_xlabel(L"$x$ (pixels)")
    ax[1].set_ylabel(L"$y$ (pixels)")
    ax[1].set_title("Voronoi Bin Map")

    # Compute the radial distance from the bin with the highest S/N
    best_bin = argmax(SN)
    r = hypot.(x̄ .- x̄[best_bin], ȳ .- ȳ[best_bin])
    ax[2].plot(hypot.(x .- x̄[best_bin], y .- ȳ[best_bin]), S./N, "k.", label=L"Input $S/N$")
    single = isone.(area)
    if any(single)
        ax[2].plot(r[single], SN[single], "x", label="Not binned")
    end
    ax[2].plot(r[.~single], SN[.~single], "o", label="Voronoi bins")
    ax[2].set_xlabel(L"$R$ (from brightest bin)")
    ax[2].set_ylabel(L"$S/N$")
    ax[2].set_xlim(extrema(r)...)
    ax[2].set_ylim(0., maximum(SN)*1.05)
    ax[2].plot(extrema(r), [target_SN, target_SN], "k--", alpha=0.5, label=L"Target $S/N$")
    ax[2].legend()

    fig, ax

end