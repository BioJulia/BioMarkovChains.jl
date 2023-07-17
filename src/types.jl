abstract type AbstractBioMarkovChain end

const LongNucOrView{N} = Union{LongSequence{<:NucleicAcidAlphabet{N}},LongSubSeq{<:NucleicAcidAlphabet{N}}}

const NUCLEICINDEXES = Dict(DNA_A => 1, DNA_C => 2, DNA_G => 3, DNA_T => 4)

const DINUCLEICINDEXES = Dict(
    LongDNA{4}("AA") => [1, 1],
    LongDNA{4}("AC") => [1, 2],
    LongDNA{4}("AG") => [1, 3],
    LongDNA{4}("AT") => [1, 4],
    LongDNA{4}("CA") => [2, 1],
    LongDNA{4}("CC") => [2, 2],
    LongDNA{4}("CG") => [2, 3],
    LongDNA{4}("CT") => [2, 4],
    LongDNA{4}("GA") => [3, 1],
    LongDNA{4}("GC") => [3, 2],
    LongDNA{4}("GG") => [3, 3],
    LongDNA{4}("GT") => [3, 4],
    LongDNA{4}("TA") => [4, 1],
    LongDNA{4}("TC") => [4, 2],
    LongDNA{4}("TG") => [4, 3],
    LongDNA{4}("TT") => [4, 4],
)

const AA20 = (AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V)

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
        bmc = new{M,I,N}(tpm, inits, n)
        return bmc
    end

    function BioMarkovChain(sequence::LongNucOrView{4}, n::Int64=1)
        tpm = transition_probability_matrix(sequence, n)
        inits = initials(sequence)
        bmc = new{Matrix{Float64},Vector{Float64},Int64}(tpm, inits, n)
        return bmc
    end
end

"""
    BMC

Alias for the type `BioMarkovChain`.
"""
const BMC = BioMarkovChain