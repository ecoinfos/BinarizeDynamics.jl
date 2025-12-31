using BinarizeDynamics
using Random

# Generate synthetic data
# N=100 sequences (approx 5000 pairs), L=100 positions
# This is a small-to-medium dataset.
N = 100
L = 100
seqs1 = [randstring("ACGT", L) for _ in 1:N]
seqs2 = [randstring("ACGT", L) for _ in 1:N]

println("Benchmarking 2.5D Analysis (Bootstrap n=100)...")
println("Dataset: N=$N sequences, L=$L positions.")

# Pre-compile
diff_dynamics(seqs1[1:10], seqs2[1:10]; test=:bootstrap, n_resamples=2, method=:phi)

# Benchmark
t = @elapsed diff_dynamics(seqs1, seqs2; test=:bootstrap, n_resamples=100, method=:phi)

println("Time for 100 resamples: $(round(t, digits=4)) seconds")
println("Estimated time for 1000 resamples: $(round(t*10, digits=4)) seconds")
