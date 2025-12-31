module BinarizeDynamicsMakieExt

using BinarizeDynamics
using Makie

"""
    plot_interaction(mat::InteractionMatrix; title="Interaction Matrix")

Plot the interaction matrix as a heatmap.
Returns a Makie Figure.
"""
function BinarizeDynamics.plot_interaction(mat::InteractionMatrix; title="Interaction Matrix")
    # Extract matrix
    M = mat.values
    L = size(M, 1)
    
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], 
        title=title,
        xlabel="Position",
        ylabel="Position",
        aspect=DataAspect()
    )
    
    # Heatmap
    hm = heatmap!(ax, 1:L, 1:L, M, colormap=:viridis)
    Colorbar(fig[1, 2], hm, label="Correlatedness")
    
    return fig
end

"""
    plot_differential(mat::InteractionMatrix; title="Differential Dynamics")

Plot the differential matrix.
"""
function BinarizeDynamics.plot_differential(mat::InteractionMatrix; title="Differential Dynamics")
    # Extract matrix
    M = mat.values
    L = size(M, 1)
    
    limit = maximum(abs.(M))
    if limit == 0
        limit = 1.0 # Avoid error on all-zero matrix
    end
    
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], 
        title=title,
        xlabel="Position",
        ylabel="Position",
        aspect=DataAspect()
    )
    
    hm = heatmap!(ax, 1:L, 1:L, M, colormap=:balance, colorrange=(-limit, limit))
    Colorbar(fig[1, 2], hm, label="Difference")
    
    return fig
end

end
