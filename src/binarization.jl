
using Base.Threads

"""
    binarize(sequences::Vector{String}; 
             positions::Union{UnitRange{Int}, Vector{Int}}=1:length(sequences[1]), 
             return_type::Symbol=:typed, 
             nthreads::Int=Threads.nthreads()) -> Union{BinarizedPairs, BitMatrix}

Convert a vector of sequences into a pairwise difference binary map.

# Arguments
- `sequences`: Vector of strings (must be equal length).
- `positions`: Indices of positions to include (default: all).
- `return_type`: :typed (returns `BinarizedPairs`) or :raw (returns `BitMatrix`).
- `nthreads`: Number of threads to use.

# Returns
- `BinarizedPairs` (default) containing the binary difference matrix, position mapping, and metadata.
- `BitMatrix` (size n_pairs x n_positions) if `return_type=:raw`.
"""
function binarize(sequences::Vector{String}; 
                  positions::Union{UnitRange{Int}, Vector{Int}}=1:length(sequences[1]), 
                  return_type::Symbol=:typed, 
                  nthreads::Int=Threads.nthreads())
    
    L_total = validate_sequences(sequences)
    
    # Validate positions
    if !all(p -> 1 <= p <= L_total, positions)
        throw(ArgumentError("Positions must be within bounds 1:$L_total"))
    end
    
    L = length(positions)
    N = length(sequences)
    mapper = PairMapper(N)
    
    # Pre-allocate output
    # We produce a BitMatrix of size (total_pairs, L)
    # Strategy: Parallelize over L columns.
    
    # 1. Convert strings to Matrix{Char} for fast random access.
    char_mat = Matrix{Char}(undef, N, L)
    for (i, seq) in enumerate(sequences)
        for (k_idx, k_pos) in enumerate(positions)
            char_mat[i, k_idx] = seq[k_pos] # Access specific position
        end
    end
    
    # 2. Parallel Processing
    # Columns of the resulting matrix
    columns = Vector{BitVector}(undef, L)
    
    @threads for k in 1:L
        # Allocate column for this position
        bv = BitVector(undef, mapper.total_pairs)
        idx = 0
        for i in 1:N
            # Loop j > i
            @inbounds for j in (i+1):N
                idx += 1
                # Difference: 1 if different, 0 if same
                bv[idx] = (char_mat[i, k] != char_mat[j, k])
            end
        end
        columns[k] = bv
    end
    
    # 3. Assemble
    data = reduce(hcat, columns)
    
    if return_type == :raw
        return data
    elseif return_type == :typed
        # Ensure positions is a Vector
        pos_vec = collect(positions)
        return BinarizedPairs(data, pos_vec, mapper)
    else
        throw(ArgumentError("Unknown return_type: $return_type. Use :typed or :raw."))
    end
end
