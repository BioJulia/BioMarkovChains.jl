using Documenter
using BioMarkovChains
# using DocThemeIndigo

# indigo = DocThemeIndigo.install(Configurations)

DocMeta.setdocmeta!(BioMarkovChains, :DocTestSetup, :(using BioMarkovChains); recursive = true)


makedocs(;
    modules = [BioMarkovChains],
    authors = "Camilo García",
    repo = "https://github.com/camilogarciabotero/BioMarkovChains.jl/blob/{commit}{path}#{line}",
    sitename = "BioMarkovChains.jl",
    format = Documenter.HTML(
        mathengine=MathJax3(),
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical="https://camilogarciabotero.github.io/BioMarkovChains.jl",
        # assets=String[indigo],
    ),
    pages = [
        "Home" => "index.md",
        "Towards Markov Chains" => "biomarkovchains.md",
        "API" => "api.md"
    ],
)

deploydocs(; repo = "https://github.com/camilogarciabotero/BioMarkovChains.jl", devbranch = "main")
