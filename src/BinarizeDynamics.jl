module BinarizeDynamics

using Statistics
using LinearAlgebra
using LoopVectorization

# Export types
export PairMapper, BinarizedData, InteractionMatrix

# Export functions
export binarize, analyze_structure, diff_dynamics, validate_sequences, plot_interaction, plot_differential

# Include source files
include("types.jl")
include("utils.jl")
include("binarization.jl")
include("analysis.jl")
include("visualization.jl")

end
