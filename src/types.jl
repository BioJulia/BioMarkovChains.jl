export BioMarkovChain, BMC, AminoAcidSeqOrView

const AminoAcidSeqOrView = Union{LongSequence{AminoAcidAlphabet}, LongSubSeq{AminoAcidAlphabet}}

abstract type AbstractBioMarkovChain end

"""
    struct BioMarkovChain{S<:DataType, M<:AbstractMatrix, I<:AbstractVector, N<:Integer} <: AbstractBioMarkovChain

A BioMarkovChain represents a Markov chain used in biological sequence analysis. It contains a transition probability matrix (tpm) and an initial distribution of probabilities (inits) and also the order of the Markov chain.

## Fields
- `alphabet::A`: The state space of the sequence whether DNA, RNA AminoAcid `DataType`s.
- `tpm::M`: The transition probability matrix.
- `inits::I`: The initial distribution of probabilities.
- `n::N`: The order of the Markov chain.

## Constructors
- `BioMarkovChain(tpm::M, inits::I, n::N=1) where {M<:AbstractMatrix, I<:AbstractVector, N<:Integer}`: Constructs a BioMarkovChain object with the provided transition probability matrix, initial distribution, and order.
- `BioMarkovChain(sequence::LongNucOrView{4}, n::Int64=1)`: Constructs a BioMarkovChain object based on the DNA sequence and transition order.

## Example

```julia
sequence = LongDNA{4}("ACTACATCTA")

model = BioMarkovChain(sequence, 2)

BioMarkovChain of DNAAlphabet{4}() and order 2:
  - Transition Probability Matrix -> Matrix{Float64}(4 × 4):
   0.4444  0.1111  0.0     0.4444
   0.4444  0.4444  0.0     0.1111
   0.0     0.0     0.0     0.0
   0.1111  0.4444  0.0     0.4444
  - Initial Probabilities -> Vector{Float64}(4 × 1):
   0.3333  0.3333  0.0     0.3333
```
"""
struct BioMarkovChain{A<:Alphabet} <: AbstractBioMarkovChain
  tpm::Matrix{Float64} # The probabilities of the TransitionProbabilityMatrix struct
  inits::Vector{Float64} # the initials distribution of probabilities
  n::Int # The order of the Markov chain
end

function BioMarkovChain(
  seq::SeqOrView{A},
  n::Int64=1
) where {A<:Alphabet}
  tpm = transition_probability_matrix(seq)
  inits = initials(seq)
  return BioMarkovChain{A}(n > 1 ? tpm^n : tpm, inits, n)
end

function BioMarkovChain(
  ::Type{A},
  tpm::Matrix{Float64},
  inits::Vector{Float64},
  n::Int64=1
) where {A<:Alphabet}
  return BioMarkovChain{A}(n > 1 ? tpm^n : tpm, inits, n)
end

typeof(i::BioMarkovChain{A}) where {A} = A

"""
    BMC

Alias for the type `BioMarkovChain`.
"""
const BMC = BioMarkovChain