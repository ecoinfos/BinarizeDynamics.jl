# BinarizeDynamics.jl

**High-performance analysis of structural dynamics in sequence data via pairwise binarization and bitwise operations.**

## Overview

BinarizeDynamics.jl abstracts sequence relationships into binary (0/1) signals using efficient `BitArray` storage, enabling high-speed 2D structural analysis and differential dynamics. It addresses the memory explosion and speed bottlenecks of standard character-based comparison methods.

## Features

- **Memory Efficiency**: ~8x reduction using `BitArray` storage.
- **High Performance**: Multithreaded binarization and bitwise analysis using `LinearAlgebra`.
- **Structural Insight**: Calculate Phi coefficients and Average Product Correction (APC) to reveal coupled mutations.
- **Visualization**: Built-in Makie recipes for interaction heatmaps.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/YourRepo/BinarizeDynamics.jl")
```

## Quick Start

```julia
using BinarizeDynamics
using Makie

# 1. Load Sequences
seqs = ["AAA", "AAB", "ABB", "BBB"]

# 2. Binarize
data = binarize(seqs)

# 3. Analyze Structure
# Computes Phi coefficient with APC correction
inter_mat = analyze_structure(data, method=:phi_apc)

# 4. Visualize
fig = plot_interaction(inter_mat; title="Topological Linkage")
display(fig)
```

## API

- `binarize(sequences)`: Convert sequences to 3D binary difference map.
- `analyze_structure(data; method)`: Compute interaction matrix.
- `diff_dynamics(mat_b, mat_a)`: Compute difference between two states.
- `plot_interaction(mat)`: Visualize the matrix.

## License

MIT License.
