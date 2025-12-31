using Documenter
using BinarizeDynamics

makedocs(
    sitename = "BinarizeDynamics.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [BinarizeDynamics],
    pages = [
        "Home" => "index.md",
        "Scientific Concepts" => "concepts.md",
        "Mathematical Background" => "math.md",
        "API Reference" => "api.md"
    ],
    remotes = nothing
)

deploydocs(
    repo = "github.com/YourRepo/BinarizeDynamics.jl.git",
)
