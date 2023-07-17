# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog and this project adheres to Semantic Versioning.

## [UNRELEASED](https://github.com/camilogarciabotero/GeneFinder.jl/compare/v0.0.10...main)

## [0.3.0]

- Refactoring the `TransitionModel` to a `BioMarkovChain` construct and better interface

## [0.2.0]

- Exapanding the alphabets to RNA and AA to its Markov chains.
- Improving runtime and allocations
- Refacotring and renaming of internal functions

## [0.1.X]

- This version has some clean ups and a new faster version of `dinucleotides` counts thanks to @gdalle. See [here the discussion](https://discourse.julialang.org/t/optimizing-dinucleotides-count-in-a-dna-sequence-type-longdna/101583/4?u=camilogarciabotero)
## [0.1.0]

- First release of `BioMarkovChains`
