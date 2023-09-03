# VoronoiBinning

[![Build Status](https://github.com/Michael-Reefe/VoronoiBinning.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Michael-Reefe/VoronoiBinning.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Introduction

This Julia package allows one to bin 2-dimensional data such that each bin has an equal (or nearly equal) signal-to-noise ratio using
Voronoi tessellations.  It is based off the python package [Vorbin](https://pypi.org/project/vorbin/) and implements most of the same
functionality, but with some minor differences (noted below). Please see [Cappellari & Copin (2003)](https://ui.adsabs.harvard.edu/abs/2003MNRAS.342..345C/abstract)
and [Diehl & Statler (2006)](https://ui.adsabs.harvard.edu/abs/2006MNRAS.368..497D/abstract) for more information about the the theory
behind the Voronoi binning and the different weighting schemes that are employed.

## Installation

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/Michael-Reefe/VoronoiBinning.jl")
```
or alternatively,
```julia
julia> ]
(Environment)> add "https://github.com/Michael-Reefe/VoronoiBinning.jl"
```

## Requirements

The package has the following dependencies:

- [LaTeXStrings](https://github.com/JuliaStrings/LaTeXStrings.jl)
- [LoopVectorization](https://github.com/JuliaSIMD/LoopVectorization.jl)
- [NearestNeighbors](https://github.com/KristofferC/NearestNeighbors.jl)
- [Printf](https://docs.julialang.org/en/v1/stdlib/Printf/) (standard library)
- [Statistics](https://docs.julialang.org/en/v1/stdlib/Statistics/) (standard library)
- [StatsBase](https://juliastats.org/StatsBase.jl/stable/) (standard library)

There is also an *optional* dependency for [PyPlot](https://github.com/JuliaPy/PyPlot.jl) if one wishes to utilize the plotting capabilities,
which allow for quick visualizations of the voronoi bins and the signal-to-noise ratio of each bin.

## Usage

There is really only one main function, with the following call signature:
```julia
voronoi2Dbinning(x, y, S, N, target_SN[, SN_func, pixel_size, bin_strategy], verbose=false, dist2_thresh=(1.2pixel_size)^2, R_thresh=0.3, cvt_tolerance=0.1, plot=false)
```
Required arguments:
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

Optional arguments:
    has a sufficiently close S/N.
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

Keyword arguments:
- `verbose::Bool=false`: A flag that enables progress and performance messages.
- `dist2_thresh::T3=1.44pixel_size^2`: Threshold on the distance^2 to consider during the bin accretion step. If a new pixel's distance^2 is within
    this threshold, then it may be accreted.
- `R_thresh::T3=0.3`: Threshold on the roundness parameter to consider during the bin accretion step. If a new pixel does not increase the roundness
    above R_thresh, then it may be accreted.
- `cvt_tolerance::T2=0.1`: The tolerance level to check against during the iterative bin adjustment procedure, if the `bin_strategy` is CentroidalVoronoi
    or WeightedVoronoi. This is a tolerance on the sum of the differences between the positions of the node generators at the current step and the previous step.
- `plot::Bool=false`: A flag that, if set, will automatically plot the voronoi bin map as well as the S/N of each bin as a function of radius.

Return values:
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

## Differences from the python package

There are a few differences that are worth mentioning:

- Due to differences in the implementation of NearestNeighbors.jl's `KDTree` and scipy's `cKDTree`, the binning results may differ slightly from the
  python version, even when using identical settings. However, the scatter in the S/N between bins should agree within the percent level or better.
  For small images (less than ~100x100) the results agree exactly.
- Because Julia uses 1-based indexing, the returned `bin_numbers` will be 1 larger than the corresponding values returned by the python version.
- The keyword arguments `dist2_thresh`, `R_thresh`, and `cvt_tolerance` have default values that match the values used by `vorbin`, but this package
  allows the user to change them if one wishes. Similarly, the weighting scheme used by the `WeightedVoronoi` defaults to the one used in `vorbin`,
  but may also be changed to be any arbitrary function of the signal and the noise.
- The `weights` return value is not the same as the `scale` return value from `vorbin`.  More specifically, `weights = 1/scale^2`.
- Due to Julia's JIT compiler and more specific optimizations implemented for this package (using LoopVectorization.jl), the runtimes compared to the 
  python package are much faster. To get a general idea, testing with the default options on my laptop with a 600x600 image runs in 224 seconds, compared 
  to 1123 seconds using the python version. For a smaller 200x200 image, the julia version ran in 4 seconds compared to 20 seconds for the python version.

## Citation

If you use this code in your research, we request that you cite the original method paper: [Cappellari & Copin (2003)](https://ui.adsabs.harvard.edu/abs/2003MNRAS.342..345C/abstract). A BibTeX-formatted entry is given below for your convenience:

```bibtex
@ARTICLE{Cappellari2003,
    author = {{Cappellari}, M. and {Copin}, Y.},
    title = "{Adaptive spatial binning of integral-field spectroscopic
        data using Voronoi tessellations}",
    journal = {MNRAS},
    eprint = {astro-ph/0302262},
    year = 2003,
    volume = 342,
    pages = {345-354},
    doi = {10.1046/j.1365-8711.2003.06541.x}
}
```
