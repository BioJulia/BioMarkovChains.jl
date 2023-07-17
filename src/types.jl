abstract type BioMarkovChain end

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
    struct TransitionModel

The `TransitionModel` struct represents a transition model used in a sequence analysis. It consists of a transition probability matrix (TransitionProbabilityMatrix) and initial distribution probabilities.

# Fields

- `TransitionProbabilityMatrix::Matrix{Float64}`: The transition probability matrix, a matrix of type Float64 representing the probabilities of transitioning from one state to another.
- `initials::Matrix{Float64}`: The initial distribution probabilities, a matrix of type Float64 representing the probabilities of starting in each state.
- `n`: is the order of the transition model, or in other words the order of the resulted Markov chain.

# Constructors

- `TransitionModel(tpm::Matrix{Float64}, initials::Matrix{Float64}; n::Int64=1)`: Constructs a `TransitionModel` object with the provided transition probability matrix and initial distribution probabilities.
"""
struct TransitionModel <: BioMarkovChain
    tpm::Matrix{Float64} # The probabilities of the TransitionProbabilityMatrix struct
    initials::Vector{Float64} # the initials distribution of probabilities
    n::Int64 # The order of the Markov chain
end