
"""
    plot_interaction(mat::InteractionMatrix; kwargs...)

Plot the interaction matrix as a heatmap.
Requires `Makie` or `CairoMakie` to be loaded.
"""
function plot_interaction(args...; kwargs...)
    error("Plotting requires Makie. Please run `using Makie`, `using GLMakie`, or `using CairoMakie` to enable this functionality.")
end

"""
    plot_differential(mat::InteractionMatrix; kwargs...)

Plot the differential matrix.
Requires `Makie` or `CairoMakie` to be loaded.
"""
function plot_differential(args...; kwargs...)
    error("Plotting requires Makie. Please run `using Makie`, `using GLMakie`, or `using CairoMakie` to enable this functionality.")
end
