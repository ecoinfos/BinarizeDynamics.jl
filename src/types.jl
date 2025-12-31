"""
    PairMapper

Maps sequence indices (i, j) to a linear pair index k for N sequences.
Indices are 1-based.
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
    pair_index(mapper, i, j)

Return the linear index k for the pair (i, j).
Expected i < j. If i > j, they are swapped.
Throws an error if i == j or indices are out of bounds.
"""
function pair_index(mapper::PairMapper, i::Int, j::Int)
    if i == j
        throw(ArgumentError("Self-pairs (i==j) are not included in the mapping."))
    end
    if i > j
        i, j = j, i
    end
    if i < 1 || j > mapper.n_sequences
        throw(BoundsError(mapper, (i, j)))
    end

    # Formula for index of pair (i,j) with i < j in a strictly upper triangular traversal:
    # row 1 has (N-1) pairs
    # row 2 has (N-2) pairs
    # ...
    # row (i-1) has (N-(i-1)) pairs
    #
    # offset = sum_{k=1}^{i-1} (N-k) = (i-1)*N - sum_{k=1}^{i-1} k
    #        = (i-1)*N - (i-1)*i/2
    # index = offset + (j - i)

    N = mapper.n_sequences
    offset = (i - 1) * N - ((i - 1) * i) ÷ 2
    return offset + (j - i)
end

"""
    pair_from_index(mapper, k)

Return the pair (i, j) corresponding to linear index k.
This is the inverse of `pair_index`.
"""
function pair_from_index(mapper::PairMapper, k::Int)
    if k < 1 || k > mapper.total_pairs
        throw(BoundsError(mapper, k))
    end

    # Inverting the offset formula can be done by solving the quadratic or binary search.
    # For N not too large, we can iteratively find i. 
    # Current max pairs for N=1000 is ~500k. 
    # An O(1) closed form solution exists but requires sqrt.

    # i(N - (i+1)/2) approx k
    # Let's use the closed form derived from triangular numbers identity:
    # m = N*(N-1)/2 - k  (this is index from end)
    # triangle root stuff.

    # Simpler approach: 
    # The row i starts at index `start_i`. 
    # start_i = (i-1)*N - i*(i-1)/2 + 1
    # We can solve for i.
    # Or just use the closed form:
    # i = N - 1 - floor(sqrt(4*N*(N-1) - 8*k + 1)/2 - 0.5) ... maybe complex.

    # Let's stick to the algebraic 1/2*i^2 ... approximation or just simple algebra.
    # k_rev = total_pairs - k
    # i_inv = floor( (sqrt(1 + 8*k_rev) - 1) / 2 )
    # i = n - 1 - i_inv
    # This acts on the reversed index.

    N = mapper.n_sequences
    # Using the "triangular root" concept on the complement
    val = 2 * (mapper.total_pairs - k)
    # m*(m+1) <= val.  Find max m.
    m = floor(Int, (sqrt(1 + 4 * val) - 1) / 2)

    i = N - 1 - m

    # Now compute start of this row to find j
    offset = (i - 1) * N - ((i - 1) * i) ÷ 2
    j = i + (k - offset)

    return (i, j)
end

"""
    BinarizedPairs

Container for pairwise binarized differences.
"""
struct BinarizedPairs
    data::BitMatrix         # Size (n_pairs, L)
    positions::Vector{Int}  # Length L, original sequence positions
    mapper::PairMapper      # Maps pairs to rows of `data`
    metadata::Dict{Symbol,Any}

    function BinarizedPairs(data::BitMatrix, positions::Vector{Int}, mapper::PairMapper; metadata=Dict{Symbol,Any}())
        n_pairs, L = size(data)

        # Validation
        if n_pairs != mapper.total_pairs
            throw(ArgumentError("Data first dimension ($n_pairs) must match total pairs in mapper ($(mapper.total_pairs))"))
        end
        if L != length(positions)
            throw(ArgumentError("Data second dimension ($L) must match length of positions vector ($(length(positions)))"))
        end

        new(data, positions, mapper, metadata)
    end
end

"""
    InteractionMatrix

Result of structural analysis (e.g., Phi coefficients).
"""
struct InteractionMatrix
    values::Matrix{Float64}   # L x L
    positions::Vector{Int}    # Length L
    method::Symbol            # e.g. :phi
    apc_applied::Bool
    metadata::Dict{Symbol,Any}

    function InteractionMatrix(values::Matrix{Float64}, positions::Vector{Int}, method::Symbol; apc_applied::Bool=false, metadata=Dict{Symbol,Any}())
        L = size(values, 1)
        if size(values, 2) != L
            throw(ArgumentError("Interaction values must be a square matrix. Got $(size(values))"))
        end
        if length(positions) != L
            throw(ArgumentError("Positions vector length ($(length(positions))) must match matrix dimension ($L)"))
        end
        new(values, positions, method, apc_applied, metadata)
    end
end

"""
    DiffResult

Result of differential analysis with statistical testing.
"""
struct DiffResult
    effect::InteractionMatrix
    pvalues::Union{Matrix{Float64}, Missing}
    test::Symbol
    n_resamples::Int
    metadata::Dict{Symbol,Any}
end
