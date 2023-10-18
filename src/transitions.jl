"""
    transition_count_matrix(sequence::LongSequence{DNAAlphabet{4}})

Compute the transition count matrix (TCM) of a given DNA sequence.

# Arguments
- `sequence::LongSequence{DNAAlphabet{4}}`: a `LongSequence{DNAAlphabet{4}}` object representing the DNA sequence.

# Returns
A `Matrix` object representing the transition count matrix of the sequence.

# Example
```
seq = LongDNA{4}("AGCTAGCTAGCT")

tcm = transition_count_matrix(seq)

4×4 Matrix{Int64}:
 0  0  3  0
 0  0  0  3
 0  3  0  0
 2  0  0  0
```
"""
function transition_count_matrix(sequence::LongNucOrView{N}) where N 
    counts = reshape(count_kmers(sequence, 2), (4,4))'
    return copy(counts)
end

function transition_count_matrix(sequence::LongAminoAcidOrView)
    matrix = [(i,j) for i in AA20, j in AA20]
    trans = transitions(sequence)
    return reshape([get(trans, t, 0) for t in matrix], size(matrix))
end

# function transition_count_matrix(sequence::LongNucOrView{N}) where N
#     alphabetsymb = eltype(sequence) == DNA ? ACGT : ACGU # could eventually use `∘(sort, unique)(sequence)` to get a very specific and sorted alphabetsymbol
#     # alphabetsymb = ∘(sort, unique)(sequence)
#     matrix = [(i,j) for i in alphabetsymb, j in alphabetsymb]
#     trans = transitions(sequence)
#     return reshape([get(trans, t, 0) for t in matrix], size(matrix))
# end

@doc raw"""
    transition_probability_matrix(sequence::LongSequence{DNAAlphabet{4}}, n::Int64=1)

Compute the transition probability matrix (TPM) of a given DNA sequence. Formally it construct `` \hat{\mathscr{M}}`` where: 
```math
\mathscr{m}_{ij} = P(X_t = j \mid X_{t-1} = i) = \frac{{P(X_{t-1} = i, X_t = j)}}{{P(X_{t-1} = i)}}
```

The transition matrices of DNA and Amino-Acids are arranged sorted and in row-wise matrices:

First the DNA matrix:

```math
\mathscr{M}_{DNA} = \begin{bmatrix}
_{AA} & _{AC} & _{AG} & _{AT} \\
_{CA} & _{CC} & _{CG} & _{CT} \\
_{GA} & _{GC} & _{GG} & _{GT} \\
_{TA} & _{TC} & _{TG} & _{TT} \\
\end{bmatrix}
```

And then, the Aminoacids:

```math
\mathscr{M}_{AA} = \begin{bmatrix}
_{AA} & _{AC} & _{AD} & \dots & _{AW} \\
_{CA} & _{CC} & _{CD} & \dots & _{CW} \\
_{DA} & _{DC} & _{DD} & \dots & _{DW} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
_{WA} & _{WC} & _{WD} & \dots & _{WW} \\
\end{bmatrix}
```

# Arguments
- `sequence::LongNucOrView{4}`: a `LongNucOrView{4}` object representing the DNA sequence.
- `n::Int64=1`: The order of the Markov model. That is the `` \hat{M}^{n}``

# Keywords

- `extended_alphabet::Bool=false`: If true will pass the extended alphabet of DNA to search

# Returns
A `Matrix` object representing the transition probability matrix of the sequence.

# Example
```
seq = dna"AGCTAGCTAGCT"

tpm = transition_probability_matrix(seq)

4×4 Matrix{Float64}:
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0
 0.0  1.0  0.0  0.0
 1.0  0.0  0.0  0.0
```
"""
function transition_probability_matrix(sequence::SeqOrView{A}, n::Int64 = 1) where A
    tcm = transition_count_matrix(sequence)
    rowsums = sum(tcm, dims=2)
    freqs = tcm ./ rowsums
    
    freqs[isnan.(freqs) .| isinf.(freqs)] .= 0.0 # Handling NaN and Inf
    
    return n > 1 ? freqs^n : freqs
end

# function transition_probability_matrix(sequence::SeqOrView{A}, n::Int64 = 1) where A
#     tcm = transition_count_matrix(sequence)
#     rowsums = sum(tcm, dims = 2)
#     freqs = tcm ./ rowsums

#     freqs[isinf.(freqs)] .= 0.0
#     freqs[isnan.(freqs)] .= 0.0

#     return freqs^(n)
# end

transition_probability_matrix(bmc::BioMarkovChain) = bmc.tpm

@testitem "tpm" begin
    using BioSequences, BioMarkovChains
    seq01 = dna"CCTCCCGGACCCTGGGCTCGGGAC"
    tpm01 = transition_probability_matrix(seq01)

    @test round.(tpm01, digits = 3) == [0.0 1.0 0.0 0.0; 0.0 0.5 0.2 0.3; 0.25 0.125 0.625 0.0; 0.0 0.667 0.333 0.0]

    # Handling NaN and Inf
    seq02 = dna"CCTCCCGGCCCTGGGCTCGGGC"
    tpm02 = transition_probability_matrix(seq02)
    @test round.(tpm02, digits = 3) == [0.0 0.0 0.0 0.0; 0.0 0.5 0.2 0.3; 0.0 0.375 0.625 0.0; 0.0 0.667 0.333 0.0]
end


@doc raw"""
    initials(sequence::SeqOrView{A}) where A

Calculate the estimated initial probabilities for a Markov chain based on a given sequence.

This function takes a sequence of states and calculates the estimated initial probabilities
of each state in the sequence for a Markov chain. The initial probabilities are estimated
by counting the occurrences of each state at the beginning of the sequence and normalizing
the counts to sum up to 1.

```math
\begin{align}
\pi{i} = P(X_{i} = i),  i \in T  \\
\sum_{i=1}^{N} \pi_{i} = 1
\end{align}
```
Now using the dinucleotides counts estimating the initials would follow:

```math
\hat{\pi_{i}} = c_{i} \sum_{k} c_{k}
```

# Arguments
- `sequence::SeqOrView{A}`: The sequence of states representing the Markov chain.

# Returns
An `Vector{Flot64}` of estimated initial probabilities for each state in the sequence.
"""
function initials(sequence::SeqOrView{A}) where A ## π̂ estimates of the initial probabilies
    inits = Array{Float64, 1}(undef, 1)
    tcm = transition_count_matrix(sequence)
    inits = sum(tcm, dims = 1) ./ sum(tcm)
    return vec(inits)
end

initials(bmc::BioMarkovChain) = bmc.inits

@doc raw"""
    sequenceprobability(sequence::LongNucOrView{4}, model::BioMarkovChain)

Compute the probability of a given sequence using a transition probability matrix and the initial probabilities distributions of a `BioMarkovModel`.

```math
P(X_1 = i_1, \ldots, X_T = i_T) = \pi_{i_1}^{T-1} \prod_{t=1}^{T-1} a_{i_t, i_{t+1}}
```

# Arguments
- `sequence::LongNucOrView{4}`: The input sequence of nucleotides.
- `tm::BioMarkovChain` is the actual data structure composed of a `tpm::Matrix{Float64}` the transition probability matrix and `initials=Vector{Float64}` the initial state probabilities.

# Returns
- `probability::Float64`: The probability of the input sequence.

# Example

```
mainseq = LongDNA{4}("CCTCCCGGACCCTGGGCTCGGGAC")
   
bmc = BioMarkovChain(mainseq)

BioMarkovChain with DNA Alphabet:
  - Transition Probability Matrix -> Matrix{Float64}(4 × 4):
   0.0     1.0     0.0     0.0
   0.0     0.5     0.2     0.3
   0.25    0.125   0.625   0.0
   0.0     0.6667  0.3333  0.0
  - Initial Probabilities -> Vector{Float64}(4 × 1):
   0.087
   0.4348
   0.3478
   0.1304
  - Markov Chain Order -> Int64:
   1

newseq = LongDNA{4}("CCTG")

    4nt DNA Sequence:
    CCTG


dnaseqprobability(newseq, bmc)
    
    0.0217
```
"""
function dnaseqprobability(
    sequence::LongNucOrView{4},
    model::BioMarkovChain
)
    init = model.inits[NUCLEICINDEXES[sequence[1]]]

    probability = init

    for t in 1:length(sequence)-1
        i, j = DINUCLEICINDEXES[@view sequence[t:t+1]]
        probability *= model.tpm[i, j]
    end
    return probability
end

# findfirst(i -> i == (AA_T, AA_R), aamatrix)

@doc raw"""
    iscoding(
        sequence::LongSequence{DNAAlphabet{4}}, 
        codingmodel::BioMarkovChain, 
        noncodingmodel::BioMarkovChain,
        η::Float64 = 1e-5
        )

Check if a given DNA sequence is likely to be coding based on a log-odds ratio.
    The log-odds ratio is a statistical measure used to assess the likelihood of a sequence being coding or non-coding. It compares the probability of the sequence generated by a coding model to the probability of the sequence generated by a non-coding model. If the log-odds ratio exceeds a given threshold (`η`), the sequence is considered likely to be coding.
    It is formally described as a decision rule:

```math
S(X) = \log \left( \frac{{P_C(X_1=i_1, \ldots, X_T=i_T)}}{{P_N(X_1=i_1, \ldots, X_T=i_T)}} \right) \begin{cases} > \eta & \Rightarrow \text{{coding}} \\ < \eta & \Rightarrow \text{{noncoding}} \end{cases}
```

# Arguments
- `sequence::LongSequence{DNAAlphabet{4}}`: The DNA sequence to be evaluated.
- `codingmodel::BioMarkovChain`: The transition model for coding regions.
- `noncodingmodel::BioMarkovChain`: The transition model for non-coding regions.
- `η::Float64 = 1e-5`: The threshold value (eta) for the log-odds ratio (default: 1e-5).

# Returns
- `true` if the sequence is likely to be coding.
- `false` if the sequence is likely to be non-coding.

# Raises
- `ErrorException`: if the length of the sequence is not divisible by 3.
- `ErrorException`: if the sequence contains a premature stop codon.

# Example

```
sequence = LondDNA("ATGGCATCTAG")
codingmodel = BioMarkovChain()
noncodingmodel = BioMarkovChain()
iscoding(sequence, codingmodel, noncodingmodel)  # Returns: true
```
"""
function iscoding(
    sequence::LongNucOrView{4},
    codingmodel::BioMarkovChain,
    noncodingmodel::BioMarkovChain,
    η::Float64 = 1e-5
)
    pcoding = dnaseqprobability(sequence, codingmodel)
    pnoncoding = dnaseqprobability(sequence, noncodingmodel)

    logodds = log(pcoding / pnoncoding)

    length(sequence) % 3 == 0 || error("The sequence is not divisible by 3")

    !hasprematurestop(sequence) || error("There is a premature stop codon in the sequence")

    if logodds > η
        return true
    else
        false
    end
end

"""
    logoddsratio(sequence::LongNucOrView{4}, model::BioMarkovChain)

Calculates the log-odds ratio between the transition probability matrix of a given DNA sequence and a reference model.

# Arguments

- `sequence::LongNucOrView{4}`: A DNA sequence or view with a length of 4 nucleotides.
- `model::BioMarkovChain`: A reference BioMarkovChain model.


# Examples

```julia
sequence = LongNucOrView{4}("ACGT")
model = BioMarkovChain(...)  # Provide appropriate initialization for BioMarkovChain
result = logoddsratio(sequence, model)
```

"""
function logoddsratio(
    sequence::LongNucOrView{4},
    model::BioMarkovChain,
)
    tpm = transition_probability_matrix(sequence)

    return log2.(tpm./model.tpm) 
end

"""
logoddsratio(model1::BioMarkovChain, model2::BioMarkovChain)

Calculates the log-odds ratio between the transition probability matrices of two BioMarkovChain models.

# Arguments 

- `model1::BioMarkovChain`: The first BioMarkovChain model.
- `model2::BioMarkovChain`: The second BioMarkovChain model.

"""
function logoddsratio(
    model1::BioMarkovChain,
    model2::BioMarkovChain,
)
    return log2.(model1.tpm ./ model2.tpm) 
end

function logoddsratioscore(
    sequence::LongNucOrView{4},
    model::BioMarkovChain,
)
    tpm = transition_probability_matrix(sequence)

    return sum(log2.(tpm./model.tpm)) / length(sequence)
end

# function logoddsratioscore.(
#     sequences::Vector{T},
#     model::BioMarkovChain,
# ) where T <: LongNucOrView{4}
#     tpm = transition_probability_matrix.(sequence)

#     return sum.(log2.(tpm./model.tpm)) / length.(sequences)
# end