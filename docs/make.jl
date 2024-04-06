using Documenter
using BioMarkovChains

DocMeta.setdocmeta!(BioMarkovChains, :DocTestSetup, :(using BioMarkovChains); recursive = true)

fmt = Documenter.HTML(
    mathengine=MathJax3(),
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical="https://camilogarciabotero.github.io/BioMarkovChains.jl",
    assets=String[],
)

pgs = [
    "Home" => "index.md",
    "Towards Markov Chains" => "biomarkovchains.md",
    "API" => "api.md"
]

makedocs(;
    modules = [BioMarkovChains],
    authors = "Camilo García",
    repo = "https://github.com/camilogarciabotero/BioMarkovChains.jl/blob/{commit}{path}#{line}",
    sitename = "BioMarkovChains.jl",
    format = fmt,
    pages = pgs,
)

deploydocs(; repo = "https://github.com/camilogarciabotero/BioMarkovChains.jl", devbranch = "main")
