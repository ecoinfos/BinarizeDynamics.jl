# Tutorial: Analyzing DNA Structural Dynamics with BinarizeDynamics.jl

Welcome to the `BinarizeDynamics.jl` tutorial. This guide will walk you through analyzing the evolutionary coupling and structural dynamics of DNA sequences. We will cover data loading, pairwise binarization, structural analysis using the Phi coefficient, and differential dynamics between two populations.

This package is designed to be efficient, even for large datasets, but we will start with a manageable example of DNA sequences to demonstrate the workflow.

# Motivation
Inspired by the pairwise representation mechanisms in modern AI models like AlphaFold2's Evoformer, BinarizeDynamics.jl seeks to identify positions where entities change in coordination. Just as detecting co-evolving residues helps predict 3D protein structures, identifying correlated variations in generic sequence data allows us to estimate where dynamic changes occur and infer the latent characteristics of a population. This principle applies not only to biological sequences like DNA mutations but also to behavioral data, such as reactivity patterns in survey responses.

This package extends the efficient binary difference method originally proposed in BinarizeComp.jl. While the predecessor focused on positional heterogeneity (1st-order statistics), BinarizeDynamics.jl generalizes this to capture high-order interactions. By leveraging this framework, we anticipate facilitating research in the following areas:

Epistatic Interaction Mapping: Deciphering complex genetic networks where gene variations at multiple loci interact non-additively to determine fitness or phenotype.

Comparative Structural Dynamics: Quantifying how evolutionary constraints shift between distinct groups (e.g., drug-resistant vs. susceptible strains) using differential dynamics analysis.

Latent Pattern Recognition in Complex Systems: Uncovering hidden dependencies between variables in non-biological categorical datasets, such as identifying coupled response patterns in large-scale social surveys.

# Prerequisites
Before starting, ensure you have Julia installed. You will need to install BinarizeDynamics and a plotting backend like CairoMakie.

```julia
using Pkg
Pkg.add(url="https://github.com/ecoinfos/BinarizeDynamics.jl") # Replace with actual registry name if registered
Pkg.add("CairoMakie") # Required for visualization
```

# Step 1: Loading Your Data
BinarizeDynamics.jl supports reading sequences directly from FASTA files or CSV files. For this tutorial, we will first generate a dummy FASTA file to simulate a dataset of aligned DNA sequences.

Note: The sequences must be Multiple Sequence Alignments (MSA), meaning all sequences must have the same length.

```julia

using BinarizeDynamics
using CairoMakie

# --- generating dummy data for the tutorial ---
# In a real scenario, you would skip this block and use your own .fasta file.
fasta_content = """
>seq1
ATGCATGCAT
>seq2
ATGCATGCAC
>seq3
ATGCGTGCAT
>seq4
ATGCGTGCAC
>seq5
CCGCGTGCAT
"""
write("example_sequences.fasta", fasta_content)
# ---------------------------------------------

# Load sequences from the FASTA file
# read_fasta returns a tuple: (headers, sequences)
headers::Vector{String}, seqs::Vector{String} = read_fasta("example_sequences.fasta")

println("Loaded $(length(seqs)) sequences of length $(length(seqs[1])).")
```

# Step 2: Pairwise Binarization
The core innovation of this package is transforming character-based sequences (A, C, G, T) into a binary difference map. Instead of looking at single positions, we look at whether pairs of sequences are the same (0) or different (1) at each position.

This process is highly optimized and multi-threaded.

```julia

# Binarize the sequences
# This creates a memory-efficient BinarizedPairs object
data::BinarizedPairs = binarize(seqs)

# You can inspect the metadata
println("Total pairwise comparisons: $(data.mapper.total_pairs)")
```

# Step 3: Structural Analysis (Phi Coefficient)
Now that we have the binary map, we can analyze the structural linkage between positions. We use the Phi coefficient, which measures the correlation between mutations at two different positions.

We also apply Average Product Correction (APC). In biological data, shared ancestry (phylogeny) creates a background signal that can obscure real interactions. APC removes this global noise to reveal direct couplings.

```julia

# Analyze the structure
# method=:phi calculates the Pearson correlation of the binary differences.
# apc=true applies the Average Product Correction to remove phylogenetic bias.
interaction_matrix::InteractionMatrix = analyze_structure(data; method=:phi, apc=true)

println("Analysis complete. Matrix size: $(size(interaction_matrix.values))")
```

# Step 4: Visualization
Visualizing the interaction matrix helps identify "modules" or clusters of positions that evolve together. We provide built-in recipes for Makie.

```julia

# Plot the interaction matrix as a heatmap
# Requires 'using CairoMakie' (or GLMakie)
fig = plot_interaction(interaction_matrix; title="Structural Couplings (Phi + APC)")

# Display the figure
display(fig)

# Save the figure to a file
save("structural_coupling.png", fig)
```

Interpretation:

High values (Yellow/Bright): Strong positive coupling. When one position mutates, the other tends to mutate as well.

Low values (Purple/Dark): Independence or weak coupling.

# Step 5: Differential Dynamics (Comparing Groups)
A powerful feature of BinarizeDynamics.jl is Differential Dynamics. This allows you to compare the structural architecture of two different groups (e.g., "Control" vs. "Treatment", or "Species A" vs. "Species B") to see how the evolutionary constraints have shifted.

Let's simulate two groups of sequences:

```julia

# Create two distinct groups of sequences
# Group 1: Sequences where pos 1 and 2 are coupled
group1::Vector{String} = ["AAAAAA", "BBAAAA", "AAAAAA", "BBAAAA"] 

# Group 2: Sequences where pos 1 and 2 are independent/random
group2::Vector{String} = ["AABBBB", "BAAAAA", "ABAAAA", "BBAAAA"]

# Perform Differential Analysis
# This calculates (Group 1 - Group 2) and runs a bootstrap test for significance.
# n_resamples=1000 is recommended for publication-quality p-values.
diff_result::DiffResult = diff_dynamics(
    group1, 
    group2; 
    test=:bootstrap, 
    n_resamples=1000, 
    method=:phi, 
    apc=true
)

# Check the results
println("Differential analysis performed with $(diff_result.n_resamples) bootstrap resamples.")
```

## Visualizing the Difference
The differential plot highlights interactions that are significantly stronger in one group compared to the other.

```julia

# Plot the differential matrix
# Red/Positive: Interaction is stronger in Group 1
# Blue/Negative: Interaction is stronger in Group 2
fig_diff = plot_differential(diff_result.effect; title="Differential Dynamics (Group 1 - Group 2)")

save("differential_dynamics.png", fig_diff)
display(fig_diff)
```

# Summary
In this tutorial, you learned how to:

1. Load DNA sequences using read_fasta.
2. Binarize the data into a difference map using binarize.
3. Analyze structure using analyze_structure with Phi and APC.
4. Compare groups using diff_dynamics with bootstrap testing.

For more mathematical details on how the Phi coefficient and APC work, please refer to the Concepts and Math pages.
