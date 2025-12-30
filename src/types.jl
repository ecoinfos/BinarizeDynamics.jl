"""
    PairMapper

Maps sequence indices (i, j) to a linear pair index k.
"""
struct PairMapper
    n_sequences::Int
    total_pairs::Int

    function PairMapper(n::Int)
        total = (n * (n - 1)) ÷ 2
        new(n, total)
    end
end

"""
    BinarizedData

Container for 3D binary difference data.
"""
struct BinarizedData
    data::BitMatrix
    mapper::PairMapper
    seq_length::Int

    function BinarizedData(data::BitMatrix, mapper::PairMapper, seq_length::Int)
        if size(data, 1) != mapper.total_pairs
            throw(ArgumentError("Data first dimension $(size(data, 1)) must match total pairs $(mapper.total_pairs)"))
        end
        if size(data, 2) != seq_length
            throw(ArgumentError("Data second dimension $(size(data, 2)) must match sequence length $seq_length"))
        end
        new(data, mapper, seq_length)
    end
end

"""
    InteractionMatrix

Result of 2D structural analysis.
"""
struct InteractionMatrix
    matrix::Matrix{Float64}
    method::Symbol
end
