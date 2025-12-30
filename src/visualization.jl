using Makie

"""
    plot_interaction(mat::InteractionMatrix; title="Interaction Matrix")

Plot the interaction matrix as a heatmap.
Returns a Makie Figure.
"""
function plot_interaction(mat::InteractionMatrix; title="Interaction Matrix")
    # Extract matrix
    M = mat.matrix
    L = size(M, 1)
    
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], 
        title=title,
        xlabel="Position",
        ylabel="Position",
        aspect=DataAspect()
    )
    
    # Heatmap
    # M is (L x L)
    # Makie heatmap expects x, y, z.
    hm = heatmap!(ax, 1:L, 1:L, M, colormap=:viridis)
    Colorbar(fig[1, 2], hm, label="Correlatedness")
    
    return fig
end

"""
    plot_differential(mat::InteractionMatrix; title="Differential Dynamics")

Plot the differential matrix.
"""
function plot_differential(mat::InteractionMatrix; title="Differential Dynamics")
    # Diverging colormap for differences
    M = mat.matrix
    L = size(M, 1)
    
    limit = maximum(abs.(M))
    
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
