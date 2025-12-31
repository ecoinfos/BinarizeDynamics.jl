# BinarizeDynamics.jl

[![Build Status](https://github.com/ecoinfos/BinarizeDynamics.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/ecoinfos/BinarizeDynamics.jl/actions/workflows/ci.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/ecoinfos/BinarizeDynamics.jl/graph/badge.svg?token=token)](https://codecov.io/gh/ecoinfos/BinarizeDynamics.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ecoinfos.github.io/BinarizeDynamics.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


**High-performance analysis of structural dynamics in sequence data via pairwise binarization and bitwise operations.**

## Overview

BinarizeDynamics.jl abstracts sequence relationships into binary (0/1) signals using efficient `BitArray` storage, enabling high-speed 2D structural analysis and differential dynamics. It addresses the memory explosion and speed bottlenecks of standard character-based comparison methods.

## Features

- **Memory Efficiency**: ~8x reduction using `BitArray` storage.
- **High Performance**: Multithreaded binarization and bitwise analysis using `LinearAlgebra`.
- **Structural Insight**: Calculate Phi coefficients and Average Product Correction (APC) to reveal coupled mutations.
- **2.5D Analysis**: Compare dynamics between groups (`diff_dynamics`) with bootstrap testing.
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
# From raw strings
seqs = ["AAA", "AAB", "ABB", "BBB"]

# From FASTA
# headers, seqs = read_fasta("data.fasta")

# From CSV/Excel (exported to CSV)
# seqs = read_csv_sequences("data.csv", "sequence_column")

# 2. Binarize
data = binarize(seqs)

# 3. Analyze Structure
# Computes Phi coefficient with APC correction
# method=:phi is default, apc=true enables correction
res = analyze_structure(data; method=:phi, apc=true)

# 4. Visualize
# Requires Makie to be loaded (e.g. `using CairoMakie`)
fig = plot_interaction(res; title="Topological Linkage")
display(fig)

# 5. 2.5D Group Comparison (Differential Dynamics)
# Compare the interaction structures of two different groups
group1 = ["AAA", "AAB", "ABB", "BBB"]
group2 = ["CCC", "CCD", "CCA", "AAA"]

# Computes difference (Group1 - Group2) with bootstrap statistical testing
diff_res = diff_dynamics(group1, group2; test=:bootstrap, n_resamples=1000)

# Visualize the shift in dynamics
fig_diff = plot_differential(diff_res.effect; title="Differential Dynamics (Group1 - Group2)")
display(fig_diff)
```

## API

- `binarize(sequences)`: Convert sequences to 3D binary difference map.
- `analyze_structure(data; method)`: Compute interaction matrix.
- `diff_dynamics(seqs1, seqs2)`: Compute differential dynamics with statistics.
- `plot_interaction(mat)`: Visualize the matrix.

## License

MIT License.
