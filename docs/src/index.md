# BinarizeDynamics.jl

*High-performance analysis of structural dynamics via pairwise binarization.*

## Overview

BinarizeDynamics.jl provides tools to analyze sequence data by converting it into binary difference maps. This allows for:

*   **Memory-efficient storage**: ~8x smaller than char matrices using `BitMatrix`.
*   **Fast analysis**: bitwise operations for massive speedups.
*   **Structural insights**: Phi coefficient and APC correction.

## Installation

```julia
using Pkg
Pkg.add("BinarizeDynamics")
```

## Quick Start

```julia
using BinarizeDynamics
# Optional: Load Makie for plotting
using CairoMakie 

# 1. Load Sequences
seqs = ["AAA", "AAB", "ABB", "BBB"]

# 2. Binarize
# Returns a BinarizedPairs object
data = binarize(seqs) 

# 3. Analyze Structure
# Returns an InteractionMatrix
res = analyze_structure(data; method=:phi, apc=true)

# 4. Visualize
# Requires Makie to be loaded
fig = plot_interaction(res)
```

## Documentation Contents

```@contents
Pages = ["concepts.md", "api.md"]
```
