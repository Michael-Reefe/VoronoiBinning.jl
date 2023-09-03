module VoronoiBinning

export BinningStrategy, 
       VoronoiStrategy, 
       InitialVoronoi, 
       CentroidalVoronoi,
       WeightedVoronoi,
       voronoi2Dbinning

using Printf
using StatsBase
using Statistics
using NearestNeighbors
using LaTeXStrings
using LoopVectorization

const DO_PLOTTING = try
    using PyPlot
    plt.plot()
    plt.close()
    true
catch
    println("WARNING: PyPlot or matplotlib not found. The plotting capability will be disabled!")
    false
end

# Some constants for setting matplotlib font sizes
const SMALL::UInt8 = 12
const MED::UInt8 = 14
const BIG::UInt8 = 16

# Matplotlib styling initialization
function __init__()
    if DO_PLOTTING
        plt.rc("font", size=MED)                   # controls default text sizes
        plt.rc("axes", titlesize=MED)              # fontsize of the axes title
        plt.rc("axes", labelsize=MED)              # fontsize of the x and y labels
        plt.rc("xtick", labelsize=SMALL)           # fontsize of the x tick labels
        plt.rc("ytick", labelsize=SMALL)           # fontsize of the y tick labels
        plt.rc("legend", fontsize=SMALL)           # legend fontsize
        plt.rc("figure", titlesize=BIG)            # fontsize of the figure title
        plt.rc("text", usetex=true)                # use LaTeX for things like axis labels
        plt.rc("font", family="Times New Roman")   # use Times New Roman font
    end
end

# Helper functions
include("voronoi_binning.jl")
include("helpers.jl")
include("plotting.jl")

end
