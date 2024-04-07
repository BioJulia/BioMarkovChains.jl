using BioMarkovChains
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(BioMarkovChains, :DocTestSetup, :(using BioMarkovChains); recursive = true)

# fmt = Documenter.HTML(
#     mathengine = MathJax3(),
#     prettyurls = get(ENV, "CI", "false") == "true",
#     canonical = "https://camilogarciabotero.github.io/BioMarkovChains.jl",
#     repolink = "https://camilogarciabotero.github.io/BioMarkovChains.jl",
#     edit_link = "main"
# )

fmt = DocumenterVitepress.MarkdownVitepress(
    repo = "https://camilogarciabotero.github.io/BioMarkovChains.jl",
    mathengine = MathJax3(),
    edit_link = "main"
)

pgs = [
    "Home" => "index.md",
    "Towards Markov Chains" => "biomarkovchains.md",
    "API" => "api.md"
]

makedocs(;
    modules = [BioMarkovChains],
    authors = "Camilo García-Botero",
    repo = "https://github.com/camilogarciabotero/BioMarkovChains.jl/",
    sitename = "BioMarkovChains.jl",
    format = fmt,
    pages = pgs,
    warnonly = true
)

deploydocs(; repo = "https://github.com/camilogarciabotero/BioMarkovChains.jl", devbranch = "main")