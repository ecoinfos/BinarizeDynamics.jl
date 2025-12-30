
using Base.Threads

"""
    binarize(sequences::Vector{String}; nthreads::Int=Threads.nthreads()) -> BinarizedData

Convert a vector of sequences into a pairwise difference binary map.
Uses multi-threading to parallelize over sequence positions.

Logic:
1. Validate sequences and convert to efficient Char/Byte matrix.
2. For each position k, compute differences for all pairs (i, j).
3. Aggregate into BinarizedData.
"""
function binarize(sequences::Vector{String}; nthreads::Int=Threads.nthreads())
    L = validate_sequences(sequences)
    N = length(sequences)
    mapper = PairMapper(N)
    
    # Pre-allocate output or use thread-local storage?
    # Strategy: Parallelize over L (positions).
    # Each position produces a column (BitVector of length total_pairs).
    # Then hcat them.
    
    # 1. Convert strings to Matrix{UInt8} or Char for fast random access.
    # Assuming codeunits is safe (UTF8 compatibility: if lengths are equal in chars, 
    # and we want char comparisons, we must be careful if unicode chars vary in byte length.
    # BUT spec says "equal length L". If using complex unicode, byte indexing might mismatch?
    # Safest: collect(String) -> Vector{Vector{Char}} or Matrix{Char}.
    # Matrix{Char} is more memory efficient than Vector{Vector{Char}}.
    
    # Let's iterate and check if ASCII.
    # For speed, Matrix{UInt8} is best if ASCII.
    # We will use `collect(Char, s)` logic into a Matrix.
    
    char_mat = Matrix{Char}(undef, N, L)
    for (i, seq) in enumerate(sequences)
        # Verify length again? valid_sequences did it.
        # Efficiently copy chars.
        # 'collect' allocates, better to iterate.
        cc = 0
        for c in seq
            cc += 1
            char_mat[i, cc] = c
        end
    end
    
    # 2. Parallel Processing
    # We produce a Vector of BitVectors, then hcat.
    columns = Vector{BitVector}(undef, L)
    
    # Helper to generate the bit column for position k
    # We can't easily closure 'char_mat' inside @threads if not careful? It's read-only, safe.
    
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
    # reduce(hcat, columns) creates a BitMatrix (Pairs x L)
    # This might trigger a large allocation, but necessary.
    data = reduce(hcat, columns)
    
    return BinarizedData(data, mapper, L)
end
