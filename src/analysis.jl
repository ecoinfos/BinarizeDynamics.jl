
"""
    analyze_structure(input::Union{BinarizedPairs, BitMatrix}; 
                      method::Symbol=:phi, 
                      apc::Bool=false, 
                      return_type::Symbol=:typed) -> Union{InteractionMatrix, Matrix{Float64}}

Calculate structural linkage between positions.

# Arguments
- `input`: `BinarizedPairs` or raw `BitMatrix` (Pairs x Positions).
- `method`: Interaction metric. Currently only `:phi` is supported.
- `apc`: Apply Average Product Correction (APC) if true.
- `return_type`: :typed (default) or :raw.
"""
function analyze_structure(input::Union{BinarizedPairs, BitMatrix}; 
                           method::Symbol=:phi, 
                           apc::Bool=false, 
                           return_type::Symbol=:typed)
    
    # 1. Unpack data
    if input isa BinarizedPairs
        X = input.data
        positions = input.positions
        n_pairs, L = size(X)
    else
        X = input
        n_pairs, L = size(X)
        positions = 1:L
    end
    
    if method == :phi_apc
        method = :phi
        apc = true
    end

    if method != :phi
        throw(ArgumentError("Only method=:phi is currently supported."))
    end

    # 2. Compute Phi Matrix
    # Phi[i, j] = (n11*n00 - n10*n01) / sqrt(n1.*n0.*n.*)
    # This is equivalent to Pearson correlation of the binary columns.
    
    counts = vec(count(X, dims=1)) # Vector{Int} of length L
    means = counts ./ n_pairs
    stds = sqrt.(means .* (1.0 .- means))
    
    # Identify constant columns (variance approx 0)
    valid_cols = stds .> 1e-9
    
    # Intersection counts (n11)
    # X is (Pairs x L), so X' * X gives (L x L) intersection counts
    # BitMatrix mul is fast.
    C = Float64.(X' * X)
    
    R = Matrix{Float64}(I, L, L) # Diagonal 1.0
    
    # Parallelize outer loop?
    # Simple double loop is usually dominated by matmul (X'*X) which is already fast.
    # Calculation of correlation from counts is O(L^2).
    
    for i in 1:L
        for j in (i+1):L
            if valid_cols[i] && valid_cols[j]
                n11 = C[i, j]
                E11 = (counts[i] * counts[j]) / n_pairs
                cov = (n11 - E11) / n_pairs
                
                denom = stds[i] * stds[j]
                val = cov / denom
                
                R[i, j] = val
                R[j, i] = val
            else
                R[i, j] = 0.0
                R[j, i] = 0.0
            end
        end
    end
    
    # 3. Apply APC if requested
    if apc
        # APC(i, j) = (mean_row[i] * mean_col[j]) / mean_total
        # Usually applied to off-diagonals.
        # We compute means of the interaction matrix (excluding diagonal).
        
        # Row means (excluding diagonal 1.0)
        row_sums = sum(R, dims=2) .- diag(R)
        row_means = row_sums ./ (L - 1)
        
        overall_mean = sum(row_means) / L
        
        if abs(overall_mean) > 1e-12
            for i in 1:L
                for j in (i+1):L
                    apc_term = (row_means[i] * row_means[j]) / overall_mean
                    R[i, j] -= apc_term
                    R[j, i] -= apc_term
                end
            end
        end
    end
    
    if return_type == :raw
        return R
    elseif return_type == :typed
        return InteractionMatrix(R, collect(positions), method; apc_applied=apc)
    else
        throw(ArgumentError("Unknown return_type: $return_type"))
    end
end

"""
    diff_dynamics(mat_b::InteractionMatrix, mat_a::InteractionMatrix; metric=:difference) -> InteractionMatrix

Compute differential dynamics (B - A).
Both matrices must correspond to the same positions.
"""
function diff_dynamics(mat_b::InteractionMatrix, mat_a::InteractionMatrix; metric::Symbol=:difference)
    if mat_b.positions != mat_a.positions
        throw(ArgumentError("Interaction matrices must be aligned (same positions) to calculate difference."))
    end
    
    if metric == :difference
        diff_vals = mat_b.values .- mat_a.values
    else
        throw(ArgumentError("Unsupported metric: $metric"))
    end
    
    return InteractionMatrix(diff_vals, mat_b.positions, :diff_dynamics; 
                             apc_applied=(mat_b.apc_applied || mat_a.apc_applied),
                             metadata=Dict{Symbol,Any}(:metric => metric))
end

"""
    diff_dynamics(seqs1::Vector{String}, seqs2::Vector{String};
                  test=:bootstrap,
                  n_resamples=1000,
                  metric=:difference,
                  method=:phi,
                  apc=false,
                  rng=Random.default_rng()) -> DiffResult

Run differential analysis with statistical testing via bootstrapping.
Resamples sequences with replacement within each group to estimate p-values.
"""
function diff_dynamics(seqs1::Vector{String}, seqs2::Vector{String}; 
                       test=:bootstrap, 
                       n_resamples=1000, 
                       metric=:difference, 
                       method=:phi,
                       apc=false,
                       rng=Random.default_rng())
    
    # 1. Compute Observe Effect
    data1 = binarize(seqs1; return_type=:typed)
    data2 = binarize(seqs2; positions=data1.positions, return_type=:typed)
    
    mat1 = analyze_structure(data1; method=method, apc=apc)
    mat2 = analyze_structure(data2; method=method, apc=apc)
    
    obs_diff = diff_dynamics(mat1, mat2; metric=metric)
    effect_matrix = obs_diff.values
    L = size(effect_matrix, 1)
    
    pvalues = nothing
    
    if test == :bootstrap
        # Bootstrap Logic
        # Test Null Hypothesis: Difference is zero (statistically significant shift)
        # We compute p-value by checking frequency where bootstrap distribution crosses 0 ???
        # Actually standard bootstrap p-value for "Difference != 0":
        # Check if 0 is in the confidence interval.
        # Two-sided p-value:
        # p = 2 * min(P(E* > 0), P(E* < 0)) normalized by total samples.
        # If mean is positive, we look at how many are < 0.
        # If mean is negative, we look at how many are > 0.
        
        # Thread-local accumulators to avoid race conditions
        # We need independent counters for each thread
        n_threads = Threads.nthreads()
        
        # Pre-allocate thread-local storage
        # each element is (count_less, count_greater) tuple of matrices
        thread_buffers = [(zeros(Int, L, L), zeros(Int, L, L)) for _ in 1:n_threads]
        
        Threads.nthreads() > 1 && println("Running bootstrap with $(Threads.nthreads()) threads...")

        Threads.@threads for k in 1:n_resamples
            tid = Threads.threadid()
            (t_less, t_greater) = thread_buffers[tid]
            
            # Resample sequence INDICES to avoid string copying overhead if possible? 
            # Vector{String} access is fast.
            # Using random indices is cleaner.
            
            idx1 = rand(rng, 1:length(seqs1), length(seqs1))
            idx2 = rand(rng, 1:length(seqs2), length(seqs2))
            
            s1_star = seqs1[idx1]
            s2_star = seqs2[idx2]
            
            # Compute raw matrices
            # Note: binarize uses @threads internally. Nested threading might happen.
            # Usually fine.
            d1_star = binarize(s1_star; positions=data1.positions, return_type=:raw)
            d2_star = binarize(s2_star; positions=data1.positions, return_type=:raw)
            
            # Analyze
            m1_star = analyze_structure(d1_star; method=method, apc=apc, return_type=:raw)
            m2_star = analyze_structure(d2_star; method=method, apc=apc, return_type=:raw)
            
            diff_star = m1_star .- m2_star
            
            # Update thread-local counts
            for i in 1:L*L
                if diff_star[i] < 0
                    t_less[i] += 1
                elseif diff_star[i] > 0
                    t_greater[i] += 1
                end
            end
        end
        
        # Reduce results from all threads
        count_less = sum(start[1] for start in thread_buffers)
        count_greater = sum(start[2] for start in thread_buffers)
        
        # Compute Two-sided P-values
        # p = 2 * min(count_less, count_greater) / n_resamples
        
        pmat = zeros(Float64, L, L)
        for i in 1:L*L
            min_count = min(count_less[i], count_greater[i])
            pmat[i] = (2.0 * min_count) / n_resamples
            # Cap at 1.0 (though logic implies <= 1.0)
            if pmat[i] > 1.0
                pmat[i] = 1.0
            end
        end
        pvalues = pmat
    elseif test == :none
        pvalues = missing
    else
         throw(ArgumentError("Unknown test mode: $test"))
    end
    
    return DiffResult(obs_diff, pvalues, test, n_resamples, Dict{Symbol,Any}())
end
