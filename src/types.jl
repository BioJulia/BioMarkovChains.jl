abstract type AbstractBioMarkovChain end

"""
    struct BioMarkovChain{S<:DataType, M<:AbstractMatrix, I<:AbstractVector, N<:Integer} <: AbstractBioMarkovChain

A BioMarkovChain represents a Markov chain used in biological sequence analysis. It contains a transition probability matrix (tpm) and an initial distribution of probabilities (inits) and also the order of the Markov chain.

# Fields
- `statespace::S`: Is the state space of the sequence whether DNA, RNA AminoAcid `DataType`s.
- `tpm::M`: The transition probability matrix.
- `inits::I`: The initial distribution of probabilities.
- `n::N`: The order of the Markov chain.

# Constructors
- `BioMarkovChain(tpm::M, inits::I, n::N=1) where {M<:AbstractMatrix, I<:AbstractVector, N<:Integer}`: Constructs a BioMarkovChain object with the provided transition probability matrix, initial distribution, and order.
- `BioMarkovChain(sequence::LongNucOrView{4}, n::Int64=1)`: Constructs a BioMarkovChain object based on the DNA sequence and transition order.

# Example
```
sequence = LongDNA{4}("ACTACATCTA")

model = BioMarkovChain(sequence, 2)
BioMarkovChain:
  - Transition Probability Matrix -> Matrix{Float64}(4 × 4):
    0.444    0.111	0.0	  0.444
    0.444    0.444	0.0	  0.111
    0.0      0.0	0.0	  0.0
    0.111    0.444	0.0	  0.444
  - Initial Probabilities -> Vector{Float64}(4 × 1):
    0.333
    0.333
    0.0
    0.333
  - Markov Chain Order:2
```
"""
struct BioMarkovChain{S<:DataType, M<:AbstractMatrix, I<:AbstractVector, N<:Integer} <: AbstractBioMarkovChain
  statespace::S # The sequence DataType (DNA,RNA,AminoAcid)
  tpm::M # The probabilities of the TransitionProbabilityMatrix struct
  inits::I # the initials distribution of probabilities
  n::N # The order of the Markov chain
  function BioMarkovChain(statespace::S, tpm::M, inits::I, n::N=1) where {S<:DataType, M<:AbstractMatrix, I<:AbstractVector, N<:Integer} 
    bmc = new{S,M,I,N}(statespace, n > 1 ? tpm^n : tpm, inits, n)
    return bmc
  end

  function BioMarkovChain(sequence::SeqOrView{A}, n::Int64=1) where A
    inits = initials(sequence)
    tpm = transition_probability_matrix(sequence)
    bmc = new{DataType,Matrix{Float64},Vector{Float64},Int64}(eltype(sequence), n > 1 ? tpm^n : tpm, inits, n)
    return bmc
  end
end

"""
    BMC

Alias for the type `BioMarkovChain`.
"""
const BMC = BioMarkovChain