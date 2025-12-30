using Documenter
using BinarizeDynamics

makedocs(
    sitename = "BinarizeDynamics.jl",
    modules = [BinarizeDynamics],
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/USERNAME/BinarizeDynamics.jl.git",
    devbranch = "main"
)
