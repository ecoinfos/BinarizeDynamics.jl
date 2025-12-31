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
    
    @testset "Math Verification" begin
        include("test_math.jl")
    end
    
    @testset "Bootstrap" begin
        include("test_bootstrap.jl")
    end

    @testset "IO Adapters" begin
        include("test_io.jl")
    end
end
