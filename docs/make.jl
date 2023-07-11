using Documenter
using BioMarkovChains
using DocThemeIndigo

indigo = DocThemeIndigo.install(Configurations)

makedocs(;
    modules = [BioMarkovChains],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://camilogarciabotero.github.io/BioMarkovChains.jl",
        assets=String[indigo],
    ),
    pages = [
        "Home" => "index.md",
        "Towards Markov Chains" => "biomarkovchains.md",
        "API" => "api.md"
    ],
    repo = "https://github.com/camilogarciabotero/BioMarkovChains.jl",
    sitename = "BioMarkovChains.jl",
)

deploydocs(; repo = "https://github.com/camilogarciabotero/BioMarkovChains.jl")
