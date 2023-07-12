"""
    transition_count_matrix(sequence::LongSequence{DNAAlphabet{4}})

Compute the transition count matrix (TCM) of a given DNA sequence.

# Arguments
- `sequence::LongSequence{DNAAlphabet{4}}`: a `LongSequence{DNAAlphabet{4}}` object representing the DNA sequence.

# Keywords

- `extended_alphabet::Bool=false`: If true will pass the extended alphabet of DNA to search

# Returns
A `TransitionCountMatrix` object representing the transition count matrix of the sequence.

# Example
```
seq = LongDNA{4}("AGCTAGCTAGCT")

tcm = transition_count_matrix(seq)

GeneFinder.TransitionCountMatrix{Dict{DNA, Int64}, Matrix{Int64}:
   A C G T
A  0 0 3 0
C  0 0 0 3
G  0 3 0 0
T  2 0 0 0

```
 """
 function transition_count_matrix(
    sequence::LongNucOrView{4};
    extended_alphabet::Bool = false
)

    transitions = dinucleotides(sequence; extended_alphabet)
    alphabetsymbols = extended_alphabet ? collect(alphabet(DNA)) : [DNA_A, DNA_C, DNA_G, DNA_T]
    tcm = TransitionCountMatrix(alphabetsymbols)
    
    for (dinucleotide, count) in transitions
        nucleotide1, nucleotide2 = dinucleotide
        i, j = tcm.order[nucleotide1], tcm.order[nucleotide2]
        tcm.counts[i, j] = count
    end

    return tcm
end


@doc raw"""
    transition_probability_matrix(sequence::LongSequence{DNAAlphabet{4}})

Compute the transition probability matrix (TPM) of a given DNA sequence. Formally it construct `` \hat{A}`` where: 
```math
a_{ij} = P(X_t = j \mid X_{t-1} = i) = \frac{{P(X_{t-1} = i, X_t = j)}}{{P(X_{t-1} = i)}}
```

# Arguments
- `sequence::LongNucOrView{4}`: a `LongNucOrView{4}` object representing the DNA sequence.
- `n::Int64=1`: The order of the Markov model. That is the `` \hat{A}^{n}``

# Keywords

- `extended_alphabet::Bool=false`: If true will pass the extended alphabet of DNA to search

# Returns
A `TPM` object representing the transition probability matrix of the sequence.

# Example
```
seq = dna"AGCTAGCTAGCT"

tpm = transition_probability_matrix(seq)

GeneFinder.TransitionProbabilityMatrix{Dict{DNA, Int64}, Matrix{Float64}:
   A   C   G   T
A  0.0 0.0 1.0 0.0
C  0.0 0.0 0.0 1.0
G  0.0 1.0 0.0 0.0
T  1.0 0.0 0.0 0.0
```
"""
function transition_probability_matrix(
    sequence::LongNucOrView{4},
    n::Int64 = 1;
    extended_alphabet::Bool = false
)
    tcm = transition_count_matrix(sequence; extended_alphabet)
    rowsums = sum(tcm.counts, dims = 2)
    freqs = tcm.counts ./ rowsums

    freqs[isinf.(freqs)] .= 0.0
    freqs[isnan.(freqs)] .= 0.0

    return TransitionProbabilityMatrix(tcm.order, freqs^(n))
end

@testitem "tpm" begin
    using BioSequences, GeneFinder
    seq = dna"CCTCCCGGACCCTGGGCTCGGGAC"
    tpm = transition_probability_matrix(seq)

    @test round.(tpm.probabilities, digits = 3) == [0.0 1.0 0.0 0.0; 0.0 0.5 0.2 0.3; 0.25 0.125 0.625 0.0; 0.0 0.667 0.333 0.0]
end

function initial_distribution(sequence::LongNucOrView{4}) ## π̂ estimates of the initial probabilies
    # initials = Array{Float64}(undef, 1, 4)
    initials = Array{Float64, 1}(undef, 1)
    counts = transition_count_matrix(sequence).counts
    initials = sum(counts, dims = 1) ./ sum(counts)
    return vec(initials)
end

"""
    transition_model(sequence::LongNucOrView{4}, n::Int64=1)

Constructs a transition model based on the given DNA sequence and transition order.

# Arguments
- `sequence::LongNucOrView{4}`: A DNA sequence represented as a `LongNucOrView{4}` object.
- `n::Int64 (optional)`: The transition order (default: 1).

# Returns
A `TransitionModel` object representing the transition model.

---

    transition_model(tpm::Matrix{Float64}, initials::Matrix{Float64}, n::Int64=1)

Builds a transtition model based on the transition probability matrix and the initial distributions. It can also calculates higer orders of the model if `n` is changed.

# Arguments
- `tpm::Matrix{Float64}`: the transition probability matrix `TPM`
- `initials::Matrix{Float64}`: the initial distributions of the model.
- `n::Int64 (optional)`: The transition order (default: 1).

# Returns
A `TransitionProbabilityMatrix` object representing the transition probability matrix.

# Example
```
sequence = LongDNA{4}("ACTACATCTA")

model = transition_model(sequence, 2)
TransitionModel:
  - Transition Probability Matrix (Size: 4 × 4):
    0.444	0.111	0.0	0.444
    0.444	0.444	0.0	0.111
    0.0	0.0	0.0	0.0
    0.111	0.444	0.0	0.444
  - Initials (Size: 1 × 4):
    0.333	0.333	0.0	0.333
  - order: 2
```
"""
function transition_model(sequence::LongNucOrView{4}, n::Int64=1)
    tpm = transition_probability_matrix(sequence, n).probabilities
    initials = initial_distribution(sequence)
    TransitionModel(tpm, initials, n)
end

function transition_model(tpm::Matrix{Float64}, initials::Matrix{Float64}, n::Int64=1)
    TransitionModel(tpm, initials, n)
end

@doc raw"""
    sequenceprobability(sequence::LongNucOrView{4}, tpm::Matrix{Float64}, initials=Vector{Float64})

Compute the probability of a given sequence using a transition probability matrix and the initial probabilities distributions.

```math
P(X_1 = i_1, \ldots, X_T = i_T) = \pi_{i_1}^{T-1} \prod_{t=1}^{T-1} a_{i_t, i_{t+1}}
```

# Arguments
- `sequence::LongNucOrView{4}`: The input sequence of nucleotides.
- `tm::TransitionModel` is the actual data structure composed of a `tpm::Matrix{Float64}` the transition probability matrix and `initials=Vector{Float64}` the initial state probabilities.

# Returns
- `probability::Float64`: The probability of the input sequence.

# Example

```
mainseq = LongDNA{4}("CCTCCCGGACCCTGGGCTCGGGAC")

tpm = transition_probability_matrix(mainseq)
    
    4×4 Matrix{Float64}:
    0.0   1.0    0.0    0.0
    0.0   0.5    0.2    0.3
    0.25  0.125  0.625  0.0
    0.0   0.667  0.333  0.0

initials = initial_distribution(mainseq)

    1×4 Matrix{Float64}:
    0.0869565  0.434783  0.347826  0.130435
    
tm = transition_model(tpm, initials)
    - Transition Probability Matrix (Size: 4 × 4):
    0.0	1.0	0.0	0.0
    0.0	0.5	0.2	0.3
    0.25	0.125	0.625	0.0
    0.0	0.667	0.333	0.0
    - Initials (Size: 1 × 4):
    0.087	0.435	0.348	0.13
    - order: 1

newseq = LondDNA("CCTG")

    4nt DNA Sequence:
    CCTG


sequenceprobability(newseq, tm)
    
    0.0217
```
"""
function sequenceprobability(
    sequence::LongNucOrView{4},
    model::TransitionModel
)
    init = model.initials[NUCLEICINDEXES[sequence[1]]]

    probability = init

    for t in 1:length(sequence)-1
        i, j = DINUCLEICINDEXES[@view sequence[t:t+1]]
        probability *= model.tpm[i, j]
    end
    return probability
end


@doc raw"""
    iscoding(
        sequence::LongSequence{DNAAlphabet{4}}, 
        codingmodel::TransitionModel, 
        noncodingmodel::TransitionModel,
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
- `codingmodel::TransitionModel`: The transition model for coding regions.
- `noncodingmodel::TransitionModel`: The transition model for non-coding regions.
- `η::Float64 = 1e-5`: The threshold value for the log-odds ratio (default: 1e-5).

# Returns
- `true` if the sequence is likely to be coding.
- `false` if the sequence is likely to be non-coding.

# Raises
- `ErrorException`: if the length of the sequence is not divisible by 3.
- `ErrorException`: if the sequence contains a premature stop codon.

# Example

```
sequence = LondDNA("ATGGCATCTAG")
codingmodel = TransitionModel()
noncodingmodel = TransitionModel()
iscoding(sequence, codingmodel, noncodingmodel)  # Returns: true
```
"""
function iscoding(
    sequence::LongNucOrView{4},
    codingmodel::TransitionModel,
    noncodingmodel::TransitionModel,
    η::Float64 = 1e-5
)
    pcoding = sequenceprobability(sequence, codingmodel)
    pnoncoding = sequenceprobability(sequence, noncodingmodel)

    logodds = log(pcoding / pnoncoding)

    length(sequence) % 3 == 0 || error("The sequence is not divisible by 3")

    !hasprematurestop(sequence) || error("There is a premature stop codon in the sequence")

    if logodds > η
        return true
    else
        false
    end
end