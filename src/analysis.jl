
"""
    analyze_structure(data::BinarizedData; method::Symbol=:phi_apc) -> InteractionMatrix

Calculate structural linkage between positions.
"""
function analyze_structure(bdata::BinarizedData; method::Symbol=:phi_apc)::InteractionMatrix
    # data matrix is (Pairs x L)
    # Rows are observations (pairs), Cols are variables (positions).
    X = bdata.data
    n_pairs, L = size(X)
    
    # Compute correlation matrix
    # We want pairwise correlations between columns.
    # Phi coefficient is Pearson correlation for binary data.
    # cor(X) in Statistics handles this efficiently for matrix?
    # BitMatrix support in Statistics.cor might be limited or Convert to float.
    # X is BitMatrix. Statistics.cor(X) returns distinct type?
    # Let's verify or implement manually for speed/memory.
    # Converting X to Matrix{Float64} is huge (8x memory). Avoid.
    
    # We can compute cor elementwise.
    # cor(A, B) = (mean(A*B) - mean(A)mean(B)) / (std(A)std(B))
    # For binary:
    # mean(A) = p_A = count(A)/n
    # var(A) = p_A * (1 - p_A)
    # mean(A*B) = count(A & B) / n
    
    # Precompute means and inv_stds
    # col_counts = sum(X, dims=1) # Returns 1xL Int array?
    # bit 'sum' or 'count'
    
    counts = vec(count(X, dims=1)) # Vector{Int} of length L
    means = counts ./ n_pairs
    stds = sqrt.(means .* (1.0 .- means))
    
    # Handle zero variance (std=0 -> const column)
    # mask where std > 0
    valid_cols = stds .> 1e-9
    
    # Result matrix
    R = Matrix{Float64}(I, L, L) # Identity diagonal
    
    # Compute upper triangle
    # Optimized: Iterate cols?
    # Or use matrix multiplication?
    # X^T * X gives the count of intersections (A & B).
    # X is (P x L). X^T is (L x P).
    # Res = X' * X is (L x L) matrix of "intersection counts".
    # This works for BitMatrix! LinearAlgebra has optimized bitmatmul?
    # Actually, X' * X with BitMatrix inputs -> Matrix{Int}.
    # This is VERY FAST (popcount).
    
    # So:
    # 1. Intersection matrix C = X' * X
    # 2. Phi[i, j] = (C[i,j]/n - p[i]*p[j]) / (std[i]*std[j])
    
    # Note: BitMatrix multiplication matches LinearAlgebra norms.
    
    C = Matrix{Int}(undef, L, L)
    # Use LinearAlgebra mul! or plain *
    # check if X' * X is supported for BitMatrix
    # It is standard in Julia.
    C = Float64.(X' * X) 
    
    for i in 1:L
        for j in (i+1):L
            if valid_cols[i] && valid_cols[j]
                # Covariance
                n11 = C[i, j]
                E11 = (counts[i] * counts[j]) / n_pairs
                cov = (n11 - E11) / n_pairs
                
                # Correlation
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
    
    if method == :phi_apc
        # APC Correction
        # APC(i, j) = (mean_row[i] * mean_col[j]) / mean_total
        # Adjusted = Raw - APC
        # Set diagonal to ?? Usually ignored. Keep 1.0 or NaN.
        
        # Calculate background means (excluding diagonal?)
        # Standard APC excludes diagonal.
        
        row_means = sum(R, dims=2) .- 1.0 # Subtract diagonal 1.0
        row_means ./= (L - 1)
        
        overall_mean = sum(row_means) / L 
        
        # If overall_mean is 0, APC is 0
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
    
    return InteractionMatrix(R, method)
end

"""
    diff_dynamics(mat_b::InteractionMatrix, mat_a::InteractionMatrix) -> InteractionMatrix

Compute differential dynamics (B - A).
"""
function diff_dynamics(mat_b::InteractionMatrix, mat_a::InteractionMatrix)
    if size(mat_b.matrix) != size(mat_a.matrix)
        throw(DimensionMismatch("Matrices must have same dimensions"))
    end
    
    # Simple difference
    # Maybe add noise filtering later
    diff = mat_b.matrix .- mat_a.matrix
    
    return InteractionMatrix(diff, :diff)
end
