# abstract type AbstractBioMarkovChain <: AbstractDiscreteMarkovChain end
abstract type AbstractBioMarkovChain end

"""
    struct BioMarkovChain{M<:AbstractMatrix, I<:AbstractVector, N<:Integer} <: AbstractBioMarkovChain

A BioMarkovChain represents a Markov chain used in biological sequence analysis. It contains a transition probability matrix (tpm) and an initial distribution of probabilities (inits) and also the order of the Markov chain.

# Fields
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
struct BioMarkovChain{M<:AbstractMatrix, I<:AbstractVector, N<:Integer} <: AbstractBioMarkovChain
    tpm::M # The probabilities of the TransitionProbabilityMatrix struct
    inits::I # the initials distribution of probabilities
    n::N # The order of the Markov chain
    function BioMarkovChain(tpm::M, inits::I, n::N=1) where {M<:AbstractMatrix, I<:AbstractVector, N<:Integer}
        bmc = new{M,I,N}(n > 1 ? tpm^n : tpm, inits, n)
        return bmc
    end

    function BioMarkovChain(sequence::SeqOrView{A}, n::Int64=1) where A
        inits = Array{Float64, 1}(undef, 1)
        tcm = transition_count_matrix(sequence)
        inits = vec(sum(tcm, dims = 1) ./ sum(tcm))

        rowsums = sum(tcm, dims = 2)
        freqs = tcm ./ rowsums

        freqs[isnan.(freqs) .| isinf.(freqs)] .= 0.0 # Handling NaN and Inf

        bmc = new{Matrix{Float64},Vector{Float64},Int64}(n > 1 ? freqs^n : freqs, inits, n)
        return bmc
    end
end

"""
    BMC

Alias for the type `BioMarkovChain`.
"""
const BMC = BioMarkovChain