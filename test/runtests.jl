using Test
using BinarizeDynamics

@testset "BinarizeDynamics.jl" begin
    @testset "Types and Utils" begin
        include("test_types_utils.jl")
    end
    
    @testset "Binarization" begin
        include("test_binarization.jl")
    end
    @testset "Analysis" begin
        include("test_analysis.jl")
    end
end
