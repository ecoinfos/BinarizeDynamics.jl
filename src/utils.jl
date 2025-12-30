"""
    validate_sequences(sequences::Vector{String})

Validates that all sequences have equal length and the vector is not empty.
Returns the length of the sequences.
"""
function validate_sequences(sequences::Vector{String})::Int
    if isempty(sequences)
        throw(ArgumentError("Sequence vector is empty"))
    end
    
    L = length(sequences[1])
    for (i, seq) in enumerate(sequences)
        if length(seq) != L
            throw(ArgumentError("Sequence at index $i has length $(length(seq)), expected $L"))
        end
    end
    return L
end

"""
    index_to_pair(k::Int, n::Int)

Converts a linear index `k` to a pair `(i, j)` for `n` items,
where `1 <= i < j <= n`.
This is the inverse of the lexicographical mapping.
"""
function index_to_pair(k::Int, n::Int)
    # Simple linear scan for now, or mathematical inversion
    # For n=4: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
    # k=1 -> (1,2)
    # k=2 -> (1,3)
    # ...
    # This can be optimized to O(1) but O(N) is fine for typically small N lookup
    # Actually, we likely need this for visualization, so O(1) is better.
    
    # Mathematical inversion for lexicographical order:
    # k = (i-1)*n - i*(i-1)/2 + (j-i)
    # This involves solving a quadratic for i.
    
    # Use robust method:
    # cp = 0
    # for i in 1:n-1
    #   rem = n - i
    #   if k <= cp + rem
    #       j = i + (k - cp)
    #       return (i, j)
    #   end
    #   cp += rem
    # end
    
    if k < 1 || k > (n * (n - 1)) ÷ 2
        throw(BoundsError("Index $k out of bounds for $n items"))
    end
    
    # optimized approach
    # Derived from: k ≈ i*n - i^2/2
    # i ≈ n - 0.5 - sqrt((n-0.5)^2 - 2k)
    # We can use this as a guess
    
    b = 2*n - 1
    # This is an approximation, let's stick to the loop for safety unless profiled
    # For very large N, the loop is O(N), which is acceptable for 'index_to_pair' usually called in loops over pairs? 
    # Actually if iterating over all k, we shouldn't call this.
    
    offset = 0
    for i in 1:n-1
        count = n - i
        if k <= offset + count
            j = i + (k - offset)
            return (i, j)
        end
        offset += count
    end
    error("Should not be reached")
end

"""
    pair_to_index(i::Int, j::Int, n::Int)

Converts a pair `(i, j)` to a linear index `k`.
Allows `i > j` (swaps automatically) but `i == j` is invalid.
"""
function pair_to_index(i::Int, j::Int, n::Int)
    if i == j
        throw(ArgumentError("i and j cannot be equal"))
    end
    if i > j
        i, j = j, i
    end
    if i < 1 || j > n
        throw(BoundsError("Indices ($i, $j) out of bounds for $n items"))
    end
    
    # Formula: sum of pairs for rows 1 to i-1, plus offset in row i
    # Pairs in row r: (n-r)
    # Sum_{r=1}^{i-1} (n-r) = n*(i-1) - i*(i-1)/2
    
    return (i - 1) * n - (i * (i - 1)) ÷ 2 + (j - i)
end
