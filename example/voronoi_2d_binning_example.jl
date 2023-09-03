using VoronoiBinning

# Note: CSV and DataFrames are not dependencies of the VoronoiBinning.jl package itself,
#       make sure you have them installed before trying to run this file!
using CSV, DataFrames

"""
    main()

This function provides a usage example for VoronoiBinning.jl using the same example as the 
python version of VorBin. This function assumes the "voronoi_2d_binning_example_input.txt" file
exists in the same directory as this file itself, and produces an output file, 
"voronoi_2d_binning_example_output.txt", also in the same directory. This test is small enough
that the differences between NearestNeighbor.jl's KDTree and scipy's cKDTree are not enough to
cause any deviations between the results using the python VorBin package and VoronoiBinning.jl.
The results have been tested and are reproduced exactly between both versions.

"""
function main()

    # Read in the data into a DataFrame object from the text file
    data = CSV.read("voronoi_2d_binning_example_input.txt", DataFrame, delim=' ', ignorerepeated=true, comment="#", header=["X", "Y", "S", "N"])
    x = data[!, :X]
    y = data[!, :Y]
    S = data[!, :S]
    N = data[!, :N]

    # Perform the binning and plot the results
    # If you don't have PyPlot.jl, this will still work, but it just wont be plotted and `fig` and `ax` will be `nothing`
    bin_numbers, xnode, ynode, xbar, ybar, SN, area, weights, fig, ax = voronoi2Dbinning(x, y, S, N, 50.0, plot=true, verbose=true)

    # save the results to an output text file
    CSV.write("voronoi_2d_binning_example_output.txt", DataFrame(x=x, y=y, N=bin_numbers))

    return bin_numbers
end


# If this script is called directly from the command line, run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
