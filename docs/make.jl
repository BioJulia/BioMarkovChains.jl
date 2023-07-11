using Documenter
using BioMarkovChains
using DocThemeIndigo

indigo = DocThemeIndigo.install(Configurations)

makedocs(;
    modules = [BioMarkovChains],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://Camilo García.github.io/BioMarkovChains.jl",
        assets=String[indigo],
    ),
    pages = [
        "Home" => "index.md",
    ],
    repo = "https://github.com/Camilo García/BioMarkovChains.jl",
    sitename = "BioMarkovChains.jl",
)

deploydocs(; repo = "https://github.com/Camilo García/BioMarkovChains.jl")
