using Test
using BinarizeDynamics

@testset "PairMapper" begin
    pm = PairMapper(4)
    @test pm.n_sequences == 4
    # Pairs for 4: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4) -> Total 6
    @test pm.total_pairs == 6
end

@testset "Index utils (PairMapper)" begin
    n = 4
    pm = PairMapper(n)
    # Expected mapping:
    # 1 -> (1,2)
    # 2 -> (1,3)
    # 3 -> (1,4)
    # 4 -> (2,3)
    # 5 -> (2,4)
    # 6 -> (3,4)
    
    @test pair_from_index(pm, 1) == (1, 2)
    @test pair_from_index(pm, 3) == (1, 4)
    @test pair_from_index(pm, 4) == (2, 3)
    @test pair_from_index(pm, 6) == (3, 4)
    
    @test_throws BoundsError pair_from_index(pm, 0)
    @test_throws BoundsError pair_from_index(pm, 7)
    
    @test pair_index(pm, 1, 2) == 1
    @test pair_index(pm, 1, 4) == 3
    @test pair_index(pm, 2, 3) == 4
    @test pair_index(pm, 3, 4) == 6
    
    # Auto-swap check
    @test pair_index(pm, 2, 1) == 1
    @test_throws ArgumentError pair_index(pm, 1, 1)
end

@testset "Core Types Constructors" begin
    # BinarizedPairs
    # data: (n_pairs, L)
    # positions: Vector{Int}
    # mapper: PairMapper
    
    pm = PairMapper(3) # 3 seqs -> 3 pairs
    L = 5
    data = trues(3, L)
    pos = collect(1:L)
    
    bp = BinarizedPairs(data, pos, pm)
    @test bp.positions == pos
    @test size(bp.data) == (3, L)
    
    # Validation fail
    @test_throws ArgumentError BinarizedPairs(trues(2, L), pos, pm) # wrong pairs
    @test_throws ArgumentError BinarizedPairs(data, [1,2], pm)     # wrong pos length
    
    # InteractionMatrix
    M = zeros(L, L)
    im = InteractionMatrix(M, pos, :phi)
    @test im.method == :phi
    @test !im.apc_applied
    
    @test_throws ArgumentError InteractionMatrix(zeros(L, L+1), pos, :phi) # not square
end

@testset "Validator" begin
    good_seqs = ["ABCD", "EFGH", "IJKL"]
    @test validate_sequences(good_seqs) == 4
    
    bad_seqs = ["ABCD", "EFG"]
    @test_throws ArgumentError validate_sequences(bad_seqs)
    
    empty_seqs = String[]
    @test_throws ArgumentError validate_sequences(empty_seqs)
end
