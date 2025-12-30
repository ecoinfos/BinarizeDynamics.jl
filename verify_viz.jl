using BinarizeDynamics
using CairoMakie
using Test

# Create dummy data
seqs = ["AAA", "AAB", "ABB", "BBB"]
data = binarize(seqs)
mat = analyze_structure(data, method=:phi_apc)

# Test Plotting
fig = plot_interaction(mat; title="Test Plot")

# Save to file to verify backend works
save("test_plot.png", fig)
println("Plot saved to test_plot.png")
