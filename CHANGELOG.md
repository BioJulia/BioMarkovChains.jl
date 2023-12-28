# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog and this project adheres to Semantic Versioning.

## [UNRELEASED](https://github.com/camilogarciabotero/GeneFinder.jl/compare/v0.0.10...main)

## [0.9.0]

- `BioMarkoChain` now has a compliant `BioSequences` alphabet.

## [0.8.1]

- Fix `BioMarkoChain` checks compats.

## [0.8.0]

- `BioMarkoChain` checks whether its row-stochastic.
- New extended functions from the `DiscreteMarkovChains` packages are now available.
- `dnaseqprobability` is much faster.
- New methods to calculate the `logg-odds-ratio-matrix` and `logg-odds-ratio-score` of a `BioSequences`.
- Improvementes on the docs.

## [0.7.0]

- Improve `BioMarkoChain` struct to be more flexible and distinguish DNA, RNA, AminoAcid #17

## [0.6.0]

- Update docs version to `1.0`.
- Update docs logo.
- Add new methods for to calculate the log-odds ratio of a sequence and a model.
- Extend new methods from the booleans/predicates of `DiscreteMarkovChains.jl` package.
- Clean code and make it more lightweight by using weakdeps.

## [0.5.0]

- Improve count of transitions 10X thanks to `VectorizedKmers.jl`.

## [0.4.1]

- Internals changes, clean code and make it more flexible.

## [0.4.0]

- Clean documentation.
- Add a `randbmc` methods.
- Improve `BMC` printing.

## [0.3.0]

- Refactoring the `TransitionModel` to a `BioMarkovChain` construct and better interface.

## [0.2.0]

- Exapanding the alphabets to RNA and AA to its Markov chains.
- Improving runtime and allocations.
- Refacotring and renaming of internal functions.

## [0.1.X]

- This version has some clean ups and a new faster version of `dinucleotides` counts thanks to @gdalle. See [here the discussion](https://discourse.julialang.org/t/optimizing-dinucleotides-count-in-a-dna-sequence-type-longdna/101583/4?u=camilogarciabotero).

## [0.1.0]

- First release of `BioMarkovChains`
