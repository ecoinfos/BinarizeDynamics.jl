module BinarizeDynamics

using Statistics
using LinearAlgebra
using Random

# Export types
export PairMapper, BinarizedData, BinarizedPairs, InteractionMatrix, DiffResult

# Export functions
export binarize, analyze_structure, diff_dynamics, validate_sequences, plot_interaction, plot_differential, read_fasta, read_csv_sequences
export pair_index, pair_from_index, index_to_pair, pair_to_index

# Include source files
include("types.jl")
include("utils.jl")
include("binarization.jl")
include("analysis.jl")
include("visualization.jl")
include("io_adapters.jl")
using .IOAdapters

end
