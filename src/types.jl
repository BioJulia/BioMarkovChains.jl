
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

# const DINUCLEOTIDES = vec([LongSequence{DNAAlphabet{2}}([i, j]) for i in ACGT, j in ACGT])

# const DINUCLEOTIDES = let
#     x = [LongSequence{DNAAlphabet{2}}([i, j]) for i in ACGT, j in ACGT]
#     # x = [string(i, j) for i in ACGT, j in ACGT]
#     Dict(x .=> 0)
# end
# const DINUCLEOTIDES = vec([LongSequence{DNAAlphabet{2}}([i, j]) for i in ACGT, j in ACGT])
# const DINUCLEOTIDES_DICT = Dict(BioMarkovChains.DINUCLEOTIDES .=> 0)

#const EXTENDED_DINUCLEOTIDES = vec([LongSequence{DNAAlphabet{4}}([i, j]) for i in alphabet(DNA), j in alphabet(DNA)])
#const EXTENDED_DINUCLEOTIDES_DICT = Dict(BioMarkovChains.EXTENDED_DINUCLEOTIDES .=> 0)

# const DINUCLEOTIDES = Dict((i, j) => 0 for (i, j) in Iterators.product(ACGT, ACGT))
# const EXTENDED_DINUCLEOTIDES = Dict((i, j) => 0 for (i, j) in Iterators.product(alphabet(DNA), alphabet(DNA)))

"""
TransitionCountMatrix(alphabet::Vector{DNA})

A data structure for storing a DNA Transition Count Matrix (TCM). The TCM is a square matrix where each row and column corresponds to a nucleotide in the given `alphabet`. The value at position (i, j) in the matrix represents the number of times that nucleotide i is immediately followed by nucleotide j in a DNA sequence. 

Fields:
- `order::Dict{DNA, Int64}`: A dictionary that maps each nucleotide in the `alphabet` to its corresponding index in the matrix.
- `counts::Matrix{Int64}`: The actual matrix of counts.

Internal function:
- `TransitionCountMatrix(alphabet::Vector{DNA})`: Constructs a new `TCM` object with the given `alphabet`. This function initializes the `order` field by creating a dictionary that maps each nucleotide in the `alphabet` to its corresponding index in the matrix. It also initializes the `counts` field to a matrix of zeros with dimensions `len x len`, where `len` is the length of the `alphabet`.

Example usage:
```julia
alphabet = [DNA_A, DNA_C, DNA_G, DNA_T]
dtcm = TCM(alphabet)
```
"""
struct TransitionCountMatrix
    order::Dict{DNA,Int64}
    counts::Matrix{Int64}

    function TransitionCountMatrix(alphabet::Vector{DNA})

        len = length(alphabet)

        order = Dict{DNA,Int}()
        for (i, nucleotide) in enumerate(sort(alphabet))
            order[nucleotide] = i
        end
        counts = zeros(Int64, len, len)
        new(order, counts)
    end
end


"""
    TransitionProbabilityMatrix(alphabet::Vector{DNA})

A data structure for storing a DNA Transition Probability Matrix (TransitionProbabilityMatrix). The TransitionProbabilityMatrix is a square matrix where each row and column corresponds to a nucleotide in the given `alphabet`. The value at position (i, j) in the matrix represents the probability of transitioning from nucleotide i to nucleotide j in a DNA sequence. 

Fields:
- `order::Dict{DNA, Int64}`: A dictionary that maps each nucleotide in the `alphabet` to its corresponding index in the matrix.
- `probabilities::Matrix{Float64}`: The actual matrix of probabilities.

Example usage:
```julia
alphabet = [DNA_A, DNA_C, DNA_G, DNA_T]
dTransitionProbabilityMatrix = TransitionProbabilityMatrix(alphabet)
```
"""
struct TransitionProbabilityMatrix
    order::Dict{DNA,Int64}
    probabilities::Matrix{Float64}
end

"""
    struct TransitionModel

The `TransitionModel` struct represents a transition model used in a sequence analysis. It consists of a transition probability matrix (TransitionProbabilityMatrix) and initial distribution probabilities.

# Fields

- `TransitionProbabilityMatrix::Matrix{Float64}`: The transition probability matrix, a matrix of type Float64 representing the probabilities of transitioning from one state to another.
- `initials::Matrix{Float64}`: The initial distribution probabilities, a matrix of type Float64 representing the probabilities of starting in each state.
- `n`: is the order of the transition model, or in other words the order of the resulted Markov chain.

# Constructors

- `TransitionModel(TransitionProbabilityMatrix::Matrix{Float64}, initials::Matrix{Float64})`: Constructs a `TransitionModel` object with the provided transition probability matrix `TransitionProbabilityMatrix` and initial distribution probabilities `initials`.
- `TransitionModel(sequence::LongSequence{DNAAlphabet{4}})`: Constructs a `TransitionModel` object based on a given DNA sequence. The transition probability matrix is calculated using `transition_probability_matrix(sequence).probabilities`, and the initial distribution probabilities are calculated using `initial_distribution(sequence)`.

"""
struct TransitionModel <: BioMarkovChain
    tpm::Matrix{Float64} # The probabilities of the TransitionProbabilityMatrix struct
    initials::Vector{Float64}
    n::Int64 # The order of the Markov chain
end
