using Test
using BinarizeDynamics

@testset "Binarization" begin
    @testset "Simple Cases" begin
        # TC-01: Minimal Case
        seqs = ["A", "B"]
        # N=2 -> 1 pair (1,2)
        # L=1
        
        # Default Typed
        data = binarize(seqs)
        @test data isa BinarizedPairs
        @test data.mapper.total_pairs == 1
        @test size(data.data) == (1, 1)
        @test data.data[1, 1] == 1 # Different
        
        # Raw Return
        raw_data = binarize(seqs, return_type=:raw)
        @test raw_data isa BitMatrix
        @test raw_data == data.data

        # TC-02: Identical
        seqs_same = ["A", "A"]
        data_same = binarize(seqs_same)
        @test data_same.data[1, 1] == 0 # Same
        
        # Multi-position
        seqs_multi = ["AA", "AB"]
        # Pos 1: A vs A -> 0
        # Pos 2: A vs B -> 1
        data_multi = binarize(seqs_multi)
        @test size(data_multi.data) == (1, 2)
        @test data_multi.data[1, 1] == 0
        @test data_multi.data[1, 2] == 1
    end
    
    @testset "Multi-Sequence" begin
        # N=3
        # 1: AAA
        # 2: AAB
        # 3: ABB
        # Pairs: (1,2), (1,3), (2,3)
        # Pos 1: A,A,A -> All same -> Pairs: 0, 0, 0
        # Pos 2: A,A,B -> (1,2)=0, (1,3)=1, (2,3)=1
        # Pos 3: A,B,B -> (1,2)=1, (1,3)=1, (2,3)=0
        
        seqs = ["AAA", "AAB", "ABB"]
        data = binarize(seqs)
        mat = data.data # (3 pairs x 3 positions)
        
        # Check Pos 1
        @test all(mat[:, 1] .== 0)
        
        # Check Pos 2
        # Pairs: (1,2): A vs A = 0
        #        (1,3): A vs B = 1
        #        (2,3): A vs B = 1
        # Indices: 1->(1,2), 2->(1,3), 3->(2,3)
        @test mat[1, 2] == 0
        @test mat[2, 2] == 1
        @test mat[3, 2] == 1
        
        # Check Pos 3
        # Pairs: (1,2): A vs B = 1
        #        (1,3): A vs B = 1
        #        (2,3): B vs B = 0
        @test mat[1, 3] == 1
        @test mat[2, 3] == 1
        @test mat[3, 3] == 0
    end
end
