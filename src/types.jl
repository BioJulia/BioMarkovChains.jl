
const LongNucOrView{N} = Union{LongSequence{<:NucleicAcidAlphabet{N}},LongSubSeq{<:NucleicAcidAlphabet{N}}}

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
struct TransitionModel
    tpm::Matrix{Float64}
    initials::Matrix{Float64}
    n::Int64
end
