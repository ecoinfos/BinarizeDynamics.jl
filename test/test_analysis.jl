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
        # Pairs:
        # (1,2): A!=B, A!=B -> 1, 1
        # (1,3): A==A, A==A -> 0, 0
        # ...
        # If the pattern of differences matches exactly, correlation should be 1.0.
        
        seqs = ["AA", "BB", "AA", "BB"]
        # Pos 1 col: [1, 0, 1, 1, 0, 1] (diff pattern for 4 seqs -> 6 pairs)
        # Pos 2 col: [1, 0, 1, 1, 0, 1]
        # Identical columns -> Correlation 1.0
        
        data = binarize(seqs)
        res = analyze_structure(data, method=:phi) # Use raw phi first
        
        # Diagonal should be 1.0
        @test res.matrix[1, 1] == 1.0
        @test res.matrix[2, 2] == 1.0
        
        # Off-diagonal should be 1.0
        @test res.matrix[1, 2] ≈ 1.0
        @test res.matrix[2, 1] ≈ 1.0
    end
    
    @testset "Anti-Coupling / Random" begin
        # Simple case: 
        # Pos 1: 0 1 0 1
        # Pos 2: 1 0 1 0
        # If we can construct sequences like this?
        # binarize result depends on pairs.
        # Direct check of analyze_structure with mock BinarizedData
        
        # Mock BitMatrix
        # 4 pairs, 2 positions
        # Col 1: 1 0 1 0
        # Col 2: 0 1 0 1 (Perfectly anticorrelated?)
        # Phi coeff should be -1.0
        
        bits = BitMatrix([1 0; 0 1; 1 0; 0 1])
        # n=4. n_seqs? doesn't matter for analysis, only internal consistency of math
        # Helper: analyze_structure computes stats on columns.
        
        # But we need a valid BinarizedData object.
        # PairMapper(N). N needs to match rows?
        # 4 pairs -> N approx? 
        # N=4 -> 6 pairs. 
        # N=3 -> 3 pairs.
        # Let's just use N=99 (dummy) in mapper, as long as size matches.
        # Actually PairMapper isn't used in analyze_structure, only bdata.data.
        
        mapper = PairMapper(4) # 6 pairs.
        # Let's make 6 rows.
        bits = BitMatrix([1 0; 1 0; 1 0; 0 1; 0 1; 0 1])
        # Col 1: 3 ones, 3 zeros. Mean = 0.5
        # Col 2: 3 zeros, 3 ones. Mean = 0.5
        # Intersection (1&1): 0.
        # Phi = (0 - 0.5*0.5) / (0.5*0.5) = -0.25 / 0.25 = -1.0
        
        bdata = BinarizedData(bits, mapper, 2)
        res = analyze_structure(bdata, method=:phi) # Raw phi
        
        @test res.matrix[1, 2] ≈ -1.0
    end
    
    @testset "APC" begin
        # Verify valid run only
        seqs = ["AAA", "BBB", "CCC"]
        data = binarize(seqs)
        res = analyze_structure(data, method=:phi_apc)
        @test size(res.matrix) == (3, 3)
    end
end
