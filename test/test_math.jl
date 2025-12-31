using Test
using BinarizeDynamics
using Statistics
using LinearAlgebra
using Random

# Brute force implementation for verification
function brute_force_phi(X::BitMatrix)
    n_pairs, L = size(X)
    R = zeros(L, L)
    for i in 1:L
        for j in 1:L
            # Extract columns
            col_i = X[:, i]
            col_j = X[:, j]
            
            # Pearson Correlation
            # cor(x, y) = mean((x-mx)(y-my)) / (sx*sy)
            m_i = mean(col_i)
            m_j = mean(col_j)
            s_i = std(col_i; corrected=false) # Population std usually used in Phi? 
            # Julia's std is sample std (div n-1).
            # Phi uses sample counts, which implies specific denominator.
            # analyze_structure implementation used:
            # stds = sqrt.(means .* (1.0 .- means)) which is POPULATION std for binary (div n).
            # So we must use corrected=false for verification.
            s_j = std(col_j; corrected=false)
            
            if s_i > 1e-9 && s_j > 1e-9
                cov = mean((col_i .- m_i) .* (col_j .- m_j))
                R[i, j] = cov / (s_i * s_j)
            else
                R[i, j] = 0.0
            end
        end
    end
    # Fix diagonal to 1.0 (if variance > 0)
    for i in 1:L
        if R[i,i] > 0.99
            R[i,i] = 1.0
        end
    end
    return R
end

@testset "Math Verification (Randomized)" begin
    # Test random matrices of various sizes
    for k in 1:5
        n_pairs = rand(10:100)
        L = rand(5:20)
        
        # Generate random bit matrix
        X = BitMatrix(rand(Bool, n_pairs, L))
        
        # Our optimized implementation
        # (We don't need full BinarizedPairs object for analyze_structure if input is BitMatrix?
        #  Wait, analyze_structure now requires Union{BinarizedPairs, BitMatrix}.
        #  If we pass BitMatrix, it works.)
        
        res_opt = analyze_structure(X, method=:phi, apc=false, return_type=:raw)
        
        # Brute force
        res_brute = brute_force_phi(X)
        
        # Compare
        @test res_opt ≈ res_brute atol=1e-10
    end
end
