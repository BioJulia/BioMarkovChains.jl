export initials, 
    transition_count_matrix,
    transition_probability_matrix,
    odds_ratio_matrix,
    log_odds_ratio_matrix,
    log_odds_ratio_score,
    markovprobability

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
function transition_count_matrix(seq::SeqOrView{A}) where {A<:Alphabet}
    return count_kmers(seq, 2).values'
end

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
```julia
seq = dna"AGCTAGCTAGCT"

tpm = transition_probability_matrix(seq)

4×4 Matrix{Float64}:
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0
 0.0  1.0  0.0  0.0
 1.0  0.0  0.0  0.0
```
"""
function transition_probability_matrix(seq::SeqOrView{A}, n::Int64 = 1) where {A <: Alphabet}
    tcm = transition_count_matrix(seq)
    rowsums = sum(tcm, dims=2)
    freqs = tcm ./ rowsums
    
    freqs[isnan.(freqs) .| isinf.(freqs)] .= 0.0 # Handling NaN and Inf
    
    return n > 1 ? freqs^n : freqs
end

transition_probability_matrix(bmc::BioMarkovChain{A}) where {A} = bmc.tpm

@doc raw"""
    initials(sequence::SeqOrView{A}) where A

Calculate the estimated initial probabilities for a Markov chain based on a given sequence.

This function takes a sequence of states and calculates the estimated initial probabilities
of each state in the sequence for a Markov chain. The initial probabilities are estimated
by counting the occurrences of each state at the beginning of the sequence and normalizing
the counts to sum up to 1.

```math
\begin{align}
\pi{i} &= P(X_{i} = i),  i \in T  \\
\sum_{i=1}^{N} \pi_{i} &= 1
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
function initials(seq::SeqOrView{A}) where {A <: Alphabet} ## π̂ estimates of the initial probabilies
    inits = Array{Float64, 1}(undef, 1)
    tcm = transition_count_matrix(seq)
    inits = sum(tcm, dims = 1) ./ sum(tcm)
    return vec(inits)
end

initials(bmc::BioMarkovChain{A}) where {A} = bmc.inits

# findfirst(i -> i == (AA_T, AA_R), aamatrix)

function odds_ratio_matrix(
    seq::SeqOrView{A},
    model::BioMarkovChain;
) where {A <: Alphabet}
    @assert model.alphabet == Alphabet(seq) "Sequence and model state space are inconsistent."
    tpm = transition_probability_matrix(seq)
    return tpm ./ model.tpm
end

@doc raw"""
    log_odds_ratio_matrix(model1::BioMarkovChain, model2::BioMarkovChain)

Calculates the log-odds ratio between the transition probability matrices of two BioMarkovChain models.

```math
\beta = \log \frac{P(x|\mathscr{m}_{1})}{P(x|\mathscr{m}_{2})}
```

Where $\mathscr{m}_{1}$ and $\mathscr{m}_{2}$ are the two models transition probability matrices.

# Arguments 

- `model1::BioMarkovChain`: The first BioMarkovChain model.
- `model2::BioMarkovChain`: The second BioMarkovChain model.

"""
function log_odds_ratio_matrix(
    modela::BioMarkovChain{A},
    modelb::BioMarkovChain{A};
    b::Number = 2
) where {A <: Alphabet}
    @assert eltype(modela) == eltype(modelb) "Models state spaces are inconsistent"
    @assert round.(sum(modela.tpm, dims=2)') == [1.0 1.0 1.0 1.0] "Model 1 transition probability matrix must be row-stochastic. That is, their row sums must be equal to 1."  
    @assert round.(sum(modelb.tpm, dims=2)') == [1.0 1.0 1.0 1.0] "Model 2 transition probability matrix must be row-stochastic. That is, their row sums must be equal to 1."  

    lorm = log.(b, modela.tpm ./ modelb.tpm)
    # lorm[isnan.(lorm) .| isinf.(lorm)] .= 0.0
    return lorm
end

@doc raw"""
    log_odds_ratio_score(sequence::SeqOrView{A}; modela::BioMarkovChain, modelb::BioMarkovChain, b::Number = 2)

Compute the log odds ratio score between a given sequence and a BioMarkovChain model.

```math
S(x) = \sum_{i=1}^{L} \beta_{x_{i}x} = \sum_{i=1} \log \frac{a^{\mathscr{m}_{1}}_{i-1} x_i}{a^{\mathscr{m}_{2}}_{i-1} x_i}
```

# Arguments
- `sequence::SeqOrView{A}`: A sequence of elements of type `A`.
- `modela::BioMarkovChain`: A BioMarkovChain model.
- `modelb::BioMarkovChain`: A BioMarkovChain model.
- `b::Number = 2`: The base of the logarithm used to compute the log odds ratio.

# Returns
The log odds ratio score between the sequence and the models.

# Example
"""
function log_odds_ratio_score(
    seq::SeqOrView{A};
    modela::BioMarkovChain=ECOLICDS,
    modelb::BioMarkovChain=ECOLINOCDS,
    b::Number = 2
) where {A <: Alphabet}
    @assert eltype(modela) == eltype(seq) "Sequence and model state space are inconsistent."
    @assert eltype(modelb) == eltype(seq) "Sequence and model state space are inconsistent."
    # @assert round.(sum(model.tpm, dims=2)') == [1.0 1.0 1.0 1.0] "Model transition probability matrix must be row-stochastic. That is, their row sums must be equal to 1."  
    # @assert round.(sum(tpm, dims=2)') == [1.0 1.0 1.0 1.0] "Sequence transition probability matrix must be row-stochastic. That is, their row sums must be equal to 1."  
    # lorm[isnan.(lorm) .| isinf.(lorm)] .= 0.0
    
    score = 0.0
    lorm = log.(b, modela.tpm ./ modelb.tpm)
    @inbounds for t in 1:length(seq)-1
        score += lorm[_dna_to_int(seq[t]), _dna_to_int(seq[t+1])]
    end

    return score
end

@doc raw"""
    markovprobability(sequence::LongNucOrView{4}, model::BioMarkovChain)

Compute the probability of a given sequence using a transition probability matrix and the initial probabilities distributions of a `BioMarkovModel`.

```math
P(X_1 = i_1, \ldots, X_T = i_T) = \pi_{i_1}^{T-1} \prod_{t=1}^{T-1} a_{i_t, i_{t+1}}
```

## Arguments
- `sequence::LongNucOrView{4}`: The input sequence of nucleotides.

## Keywords
- `model::BioMarkovChain=ECOLICDS`: A given `BioMarkovChain` model.
- `logscale::Bool=false`: If true, the function will return the log2 of the probability.
- `b::Number=2`: The base of the logarithm used to compute the log odds ratio.

# Returns
- `probability::Float64`: The probability of the input sequence given the model.

# Example

```julia
seq = LongDNA{4}("CGCGCGCGCGCGCGCGCGCGCGCGCG")
   
markovprobability(seq, model=CPGPOS, logscale=true)
    -45.073409957110556

markovprobability(seq, model=CPGNEG, logscale=true)
    -74.18912168395339
```
"""
function markovprobability(
    seq::NucleicSeqOrView{A};
    model::BioMarkovChain{A}=ECOLICDS,
    logscale::Bool=true
) where {A<:NucleicAcidAlphabet}
    @assert Alphabet(model) == Alphabet(seq) "Sequence and model state space are inconsistent."
    
    init = logscale ? log2(model.inits[_dna_to_int(seq[1])]) : model.inits[_dna_to_int(seq[1])]
    probability = init

    @inbounds for t in 1:length(seq)-1
        if logscale
            logmodel = log2.(model.tpm)
            probability += logmodel[_dna_to_int(seq[t]), _dna_to_int(seq[t+1])]
        else
            probability *= model.tpm[_dna_to_int(seq[t]), _dna_to_int(seq[t+1])]
        end
    end

    return probability
end