using Documenter
using MultifractalTools

# This is needed to make the plotting extension work during the build
using GLMakie 
GLMakie.activate!(inline=true) # Ensures plots are embeddable

# Set DocMeta to point to the source code repo
DocMeta.setdocmeta!(MultifractalTools, :DocTestSetup, :(using MultifractalTools); recursive=true)

makedocs(
    sitename = "MultifractalTools.jl",
    modules = [MultifractalTools],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true", # Makes clean URLs on GitHub
        canonical = "https://github.com/NicolauUP/MultifractalTools.jl", # <-- CHANGE THIS
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ]
)

# This deploys the documentation to GitHub Pages (optional, for later)
deploydocs(
    repo = "https://github.com/NicolauUP/MultifractalTools.jl", # <-- CHANGE THIS
    devbranch = "main"
)