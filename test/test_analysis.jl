using Test
using BinarizeDynamics
using LinearAlgebra

@testset "Analysis" begin
    @testset "Perfect Coupling" begin
        # Construct sequences where Pos 1 and Pos 2 are perfectly coupled.
        # e.g. whenever Pos 1 flips, Pos 2 flips.
        # Seq 1: A A
        # Seq 2: B B
        # Seq 3: A A
        # Seq 4: B B
        
        seqs = ["AA", "BB", "AA", "BB"]
        
        data = binarize(seqs)
        # New API: method=:phi, apc=false
        res = analyze_structure(data, method=:phi, apc=false)
        
        # Diagonal should be 1.0
        @test res.values[1, 1] == 1.0
        @test res.values[2, 2] == 1.0
        
        # Off-diagonal should be 1.0 (perfectly correlated)
        @test res.values[1, 2] ≈ 1.0
        @test res.values[2, 1] ≈ 1.0
    end
    
    @testset "Anti-Coupling / Random" begin
        # Mock BitMatrix manual construction
        # 6 pairs (N=4)
        mapper = PairMapper(4) 
        
        # Make 6 rows.
        # Col 1: 3 ones, 3 zeros. Mean = 0.5
        # Col 2: 3 zeros, 3 ones. Mean = 0.5
        # Intersection (1&1): 0.
        # Phi = (0 - 0.5*0.5) / (0.5*0.5) = -0.25 / 0.25 = -1.0
        
        bits = BitMatrix([1 0; 1 0; 1 0; 0 1; 0 1; 0 1])
        positions = [1, 2]
        
        # BinarizedPairs(data, positions, mapper)
        bdata = BinarizedPairs(bits, positions, mapper)
        
        res = analyze_structure(bdata, method=:phi, apc=false)
        
        @test res.values[1, 2] ≈ -1.0
    end
    
    @testset "APC" begin
        # Verify valid run only
        seqs = ["AAA", "BBB", "CCC"]
        data = binarize(seqs)
        
        # Raw Phi
        res_phi = analyze_structure(data, method=:phi, apc=false)
        
        # Phi + APC
        res_apc = analyze_structure(data, method=:phi, apc=true)
        
        # Legacy :phi_apc support
        res_legacy = analyze_structure(data, method=:phi_apc)
        
        @test size(res_apc.values) == (3, 3)
        @test res_apc.apc_applied == true
        @test res_legacy.values == res_apc.values
        
        # APC should generally reduce the average background signal.
        # But for strictly correlated seqs like this, values might shift.
        # Just check it runs and produces different values.
        @test res_phi.values != res_apc.values
    end
    
    @testset "Differential Dynamics" begin
        M1 = InteractionMatrix(ones(2,2), [1,2], :phi)
        M2 = InteractionMatrix(fill(0.5, 2, 2), [1,2], :phi)
        
        diff = diff_dynamics(M1, M2) # M1 - M2
        @test all(diff.values .== 0.5)
        @test diff.method == :diff_dynamics
    end
end
