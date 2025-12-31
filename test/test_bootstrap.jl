using Test
using BinarizeDynamics
using Random

@testset "Differential Dynamics & Bootstrap" begin
    # Create two groups of sequences
    # Group 1: High coupling between (1,2)
    # Group 2: No coupling
    
    # 6 sequences each
    group1 = ["AA", "BB", "AA", "BB", "AA", "BB"] # Perfect coupling (1,2)
    group2 = ["AB", "BA", "AA", "BB", "AB", "BA"] # Mixed, lower coupling
    
    # 1. No Test
    res_none = diff_dynamics(group1, group2; test=:none, method=:phi)
    @test res_none isa DiffResult
    @test res_none.test == :none
    @test ismissing(res_none.pvalues)
    
    # Check effect size
    # M1(1,2) should be 1.0
    # M2(1,2) should be lower.
    eff = res_none.effect.values[1, 2]
    @test eff > 0.0 # Group 1 > Group 2
    
    # 2. Bootstrap Test
    # Use small n_resamples for speed
    res_boot = diff_dynamics(group1, group2; test=:bootstrap, n_resamples=50, method=:phi, rng=Random.default_rng())
    
    @test res_boot isa DiffResult
    @test res_boot.test == :bootstrap
    @test res_boot.n_resamples == 50
    @test !ismissing(res_boot.pvalues)
    @test size(res_boot.pvalues) == (2, 2)
    
    # With clear signal, p-value should be relatively low (or at least calculated)
    p_val = res_boot.pvalues[1, 2]
    @test 0.0 <= p_val <= 1.0
    
    # 3. Validation
    @test_throws ArgumentError diff_dynamics(group1, group2; test=:invalid)
end
