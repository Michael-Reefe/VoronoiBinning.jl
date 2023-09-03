using VoronoiBinning
using Test

@testset "method tests" begin
    vorbin_default_small = voronoi2Dbinning([1,2,1,2], [1,1,2,2], [4.,2.,2.,0.5], [0.5,1.,0.5,0.25], 3.)
    @test all(voronoi2Dbinning([1,2,1,2], [1,1,2,2], [4.,2.,2.,0.5], [0.5,1.,0.5,0.25], 3., WeightedVoronoi()) .≈ 
        vorbin_default_small)
    @test all(voronoi2Dbinning([1,2,1,2], [1,1,2,2], [4.,2.,2.,0.5], [0.5,1.,0.5,0.25], 3., VoronoiBinning.default_SN_func) .≈ 
        vorbin_default_small)
    @test all(voronoi2Dbinning([1,2,1,2], [1,1,2,2], [4.,2.,2.,0.5], [0.5,1.,0.5,0.25], 3., 1.0) .≈ 
        vorbin_default_small)
    @test all(voronoi2Dbinning([1,2,1,2], [1,1,2,2], [4.,2.,2.,0.5], [0.5,1.,0.5,0.25], 3., 1.0, WeightedVoronoi()) .≈ 
        vorbin_default_small)
    @test all(voronoi2Dbinning([1,2,1,2], [1,1,2,2], [4.,2.,2.,0.5], [0.5,1.,0.5,0.25], 3., VoronoiBinning.default_SN_func, 1.0) .≈ 
        vorbin_default_small)
    @test all(voronoi2Dbinning([1,2,1,2], [1,1,2,2], [4.,2.,2.,0.5], [0.5,1.,0.5,0.25], 3., VoronoiBinning.default_SN_func, WeightedVoronoi()) .≈ 
        vorbin_default_small)
    @test all(voronoi2Dbinning([1,2,1,2], [1,1,2,2], [4.,2.,2.,0.5], [0.5,1.,0.5,0.25], 3., VoronoiBinning.default_SN_func, 1.0, WeightedVoronoi()) .≈ 
        vorbin_default_small)
end

