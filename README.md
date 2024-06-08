<div align="center">
  <img src="docs/src/assets/logo.svg" height="150"><br/>
  <i>Representing biological sequences as Markov chains</i><br/><br/>
</div>

<div align="center">

[![Documentation](https://img.shields.io/badge/documentation-online-blue.svg?logo=Julia&logoColor=white)](https://camilogarciabotero.github.io/BioMarkovChains.jl/dev/)
[![Latest Release](https://img.shields.io/github/release/camilogarciabotero/BioMarkovChains.jl.svg)](https://github.com/camilogarciabotero/BioMarkovChains.jl/releases/latest)
[![DOI](https://zenodo.org/badge/665161607.svg)](https://zenodo.org/badge/latestdoi/665161607)
<br/>
[![CI Workflow](https://github.com/camilogarciabotero/BioMarkovChains.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/camilogarciabotero/BioMarkovChains.jl/actions/workflows/CI.yml)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/camilogarciabotero/BioMarkovChains.jl/blob/main/LICENSE)
[![Work in Progress](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FBioMarkovChains&query=total_requests&suffix=%2Fmonth&label=Downloads)](http://juliapkgstats.com/pkg/BioMarkovChains)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

</div>

***

# BioMarkovChains

> A Julia package to represent biological sequences as Markov chains

## Installation

<p>
BioMarkovChains is a &nbsp;
    <a href="https://julialang.org">
        <img src="https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/master/images/julia.ico" width="16em">
        Julia Language
    </a>
    &nbsp; package. To install BioMarkovChains,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd>
    key in the REPL to use the package mode, then type the following command
</p>

```julia
pkg> add BioMarkovChains
```

## Creating Markov chain out of DNA sequences

An important step before developing several gene finding algorithms consist of having a Markov chain representation of the DNA. To do so, we implemented the `BioMarkovChain` method that will capture the initials and transition probabilities of a DNA sequence (`LongSequence`) and will create a dedicated object storing relevant information of a DNA Markov chain. Here an example:

Let find one ORF in a random `LongDNA` :

```julia
using BioSequences, BioMarkovChains, GeneFinder


# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"

orfseq = findorfs(seq)[3] |> sequence

21nt DNA Sequence:
ATGCGTCGAATGGCACGGTGA
```

If we translate it, we get a 7aa sequence:

```julia
translate(orfseq)

7aa Amino Acid Sequence:
MRRMAR*
```

Now supposing I do want to see how transitions are occurring in this ORF sequence, the I can use the `BioMarkovChain` method and tune it to 2nd-order Markov chain:

```julia
BioMarkovChain(orfdna, 2)

BioMarkovChain of DNAAlphabet{4}() and order 1:
  - Transition Probability Matrix -> Matrix{Float64}(4 × 4):
   0.25    0.25    0.0     0.5
   0.25    0.0     0.75    0.0
   0.25    0.25    0.25    0.25
   0.0     0.25    0.75    0.0
  - Initial Probabilities -> Vector{Float64}(4 × 1):
   0.2     0.2     0.4     0.2

```

This is  useful to later create HMMs and calculate sequence probability based on a given model, for instance we now have the *E. coli* CDS and No-CDS transition models or Markov chain implemented:

```julia
ECOLICDS

BioMarkovChain of DNAAlphabet{4}() and order 1:
  - Transition Probability Matrix -> Matrix{Float64}(4 × 4):
   0.31    0.224   0.199   0.268
   0.251   0.215   0.313   0.221
   0.236   0.308   0.249   0.207
   0.178   0.217   0.338   0.267
  - Initial Probabilities -> Vector{Float64}(4 × 1):
   0.245   0.243   0.273   0.239
```

What is then the probability of the previous random Lambda phage DNA sequence given this model?

```julia
dnaseqprobability(orfseq, ECOLICDS)

1.1061824843755975e-12
```

This is off course not very informative, but we can later use different criteria to then classify new ORFs. For a more detailed explanation see the [docs](https://camilogarciabotero.github.io/BioMarkovChains.jl/dev/biomarkovchains/)
