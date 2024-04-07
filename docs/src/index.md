# BioMarkovChains

```@raw html

<p align="center">
  <img src="logo.png" height="150"><br/>
  <i>Representing biological sequences as Markov chains</i>
</p>
```

***

## Overview

> This package aim to represent BioSequences types as Markov chains to perform different operations and predictions

## Installation

You can install `BioMarkovChains` from the julia REPL. Press `]` to enter pkg
mode, and enter the following:

    add BioMarkovChains

If you are interested in the cutting edge of the development, please
check out the master branch to try new features before release.

## Create a BioMarkovChain of DNA




``` julia
using BioSequences, BioMarkovChains

sequence = dna"CCTCCCGGACCCTGGGCTCGGGAC"

BioMarkovChain(sequence)
```

    BioMarkovChain of DNAAlphabet{4}() and order 1:
      - Transition Probability Matrix --> Matrix{Float64}(4 × 4):
       0.0     1.0     0.0     0.0
       0.0     0.5     0.2     0.3
       0.25    0.125   0.625   0.0
       0.0     0.6667  0.3333  0.0
      - Initial Probabilities -> Vector{Float64}(4 × 1):
       0.087
       0.4348
       0.3478
       0.1304

Note that, sometimes the dinucleotides transition do not harbor important biological meaning, whereas trinucleotides or codons are, in
fact, the building block of proteins. Therefore, sometimes the transition model we want to build is usually a second-order Markov chain, that represents the possible transitions of a trinucleotide.

A very nice nice property of the transition probability matrix is that the *n-step transition probability matrix* ``\mathscr{M}^{n} = (\mathscr{m}_{ij}(n))``, that is the *n*th power of ``\mathscr{M}`` represents ``i \rightarrow j`` transitions in *n* steps. We can also have higher order transition models as:

``` julia
BioMarkovChain(sequence, 2)
```

    BioMarkovChain of DNAAlphabet{4}() and order 2:
      - Transition Probability Matrix --> Matrix{Float64}(4 × 4):
       0.0     0.5     0.2     0.3
       0.05    0.475   0.325   0.15
       0.1562  0.3906  0.4156  0.0375
       0.0833  0.375   0.3417  0.2
      - Initial Probabilities -> Vector{Float64}(4 × 1):
       0.087
       0.4348
       0.3478
       0.1304