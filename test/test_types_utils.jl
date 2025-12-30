using Test
using BinarizeDynamics

@testset "PairMapper" begin
    pm = PairMapper(4)
    @test pm.n_sequences == 4
    # Pairs for 4: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4) -> Total 6
    @test pm.total_pairs == 6
end

@testset "Index utils" begin
    n = 4
    # Expected mapping:
    # 1 -> (1,2)
    # 2 -> (1,3)
    # 3 -> (1,4)
    # 4 -> (2,3)
    # 5 -> (2,4)
    # 6 -> (3,4)
    
    @test BinarizeDynamics.index_to_pair(1, n) == (1, 2)
    @test BinarizeDynamics.index_to_pair(3, n) == (1, 4)
    @test BinarizeDynamics.index_to_pair(4, n) == (2, 3)
    @test BinarizeDynamics.index_to_pair(6, n) == (3, 4)
    
    @test_throws BoundsError BinarizeDynamics.index_to_pair(0, n)
    @test_throws BoundsError BinarizeDynamics.index_to_pair(7, n)
    
    @test BinarizeDynamics.pair_to_index(1, 2, n) == 1
    @test BinarizeDynamics.pair_to_index(1, 4, n) == 3
    @test BinarizeDynamics.pair_to_index(2, 3, n) == 4
    @test BinarizeDynamics.pair_to_index(3, 4, n) == 6
    
    # Auto-swap check
    @test BinarizeDynamics.pair_to_index(2, 1, n) == 1
    @test_throws ArgumentError BinarizeDynamics.pair_to_index(1, 1, n)
end

@testset "Validator" begin
    good_seqs = ["ABCD", "EFGH", "IJKL"]
    @test validate_sequences(good_seqs) == 4
    
    bad_seqs = ["ABCD", "EFG"]
    @test_throws ArgumentError validate_sequences(bad_seqs)
    
    empty_seqs = String[]
    @test_throws ArgumentError validate_sequences(empty_seqs)
end
