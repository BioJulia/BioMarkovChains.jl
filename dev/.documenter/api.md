

- [`BioMarkovChains.BMC`](#BioMarkovChains.BMC)
- [`BioMarkovChains.BioMarkovChain`](#BioMarkovChains.BioMarkovChain)
- [`BioMarkovChains.dnaseqprobability`](#BioMarkovChains.dnaseqprobability-Union{Tuple{A},%20Tuple{Union{BioSequences.LongSequence{A},%20BioSequences.LongSubSeq{A}},%20BioMarkovChain}}%20where%20A)
- [`BioMarkovChains.initials`](#BioMarkovChains.initials-Union{Tuple{Union{BioSequences.LongSequence{A},%20BioSequences.LongSubSeq{A}}},%20Tuple{A}}%20where%20A)
- [`BioMarkovChains.log_odds_ratio_matrix`](#BioMarkovChains.log_odds_ratio_matrix-Union{Tuple{A},%20Tuple{Union{BioSequences.LongSequence{A},%20BioSequences.LongSubSeq{A}},%20BioMarkovChain}}%20where%20A)
- [`BioMarkovChains.log_odds_ratio_matrix`](#BioMarkovChains.log_odds_ratio_matrix-Tuple{BioMarkovChain,%20BioMarkovChain})
- [`BioMarkovChains.log_odds_ratio_score`](#BioMarkovChains.log_odds_ratio_score-Union{Tuple{A},%20Tuple{Union{BioSequences.LongSequence{A},%20BioSequences.LongSubSeq{A}},%20BioMarkovChain}}%20where%20A)
- [`BioMarkovChains.perronfrobenius`](#BioMarkovChains.perronfrobenius-Union{Tuple{Union{BioSequences.LongSequence{A},%20BioSequences.LongSubSeq{A}}},%20Tuple{A}}%20where%20A)
- [`BioMarkovChains.transition_count_matrix`](#BioMarkovChains.transition_count_matrix-Union{Tuple{Union{BioSequences.LongSequence{A},%20BioSequences.LongSubSeq{A}}},%20Tuple{A}}%20where%20A)
- [`BioMarkovChains.transition_probability_matrix`](#BioMarkovChains.transition_probability_matrix-Union{Tuple{Union{BioSequences.LongSequence{A},%20BioSequences.LongSubSeq{A}}},%20Tuple{A},%20Tuple{Union{BioSequences.LongSequence{A},%20BioSequences.LongSubSeq{A}},%20Int64}}%20where%20A)

<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='BioMarkovChains.BMC' href='#BioMarkovChains.BMC'>#</a>&nbsp;<b><u>BioMarkovChains.BMC</u></b> &mdash; <i>Type</i>.




```julia
BMC
```


Alias for the type `BioMarkovChain`.


[source](https://github.com/camilogarciabotero/BioMarkovChains.jl/)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='BioMarkovChains.BioMarkovChain' href='#BioMarkovChains.BioMarkovChain'>#</a>&nbsp;<b><u>BioMarkovChains.BioMarkovChain</u></b> &mdash; <i>Type</i>.




```julia
struct BioMarkovChain{S<:DataType, M<:AbstractMatrix, I<:AbstractVector, N<:Integer} <: AbstractBioMarkovChain
```


A BioMarkovChain represents a Markov chain used in biological sequence analysis. It contains a transition probability matrix (tpm) and an initial distribution of probabilities (inits) and also the order of the Markov chain.

**Fields**
- `alphabet::A`: Is the state space of the sequence whether DNA, RNA AminoAcid `DataType`s.
  
- `tpm::M`: The transition probability matrix.
  
- `inits::I`: The initial distribution of probabilities.
  
- `n::N`: The order of the Markov chain.
  

**Constructors**
- `BioMarkovChain(tpm::M, inits::I, n::N=1) where {M<:AbstractMatrix, I<:AbstractVector, N<:Integer}`: Constructs a BioMarkovChain object with the provided transition probability matrix, initial distribution, and order.
  
- `BioMarkovChain(sequence::LongNucOrView{4}, n::Int64=1)`: Constructs a BioMarkovChain object based on the DNA sequence and transition order.
  

**Example**

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



[source](https://github.com/camilogarciabotero/BioMarkovChains.jl/)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='BioMarkovChains.dnaseqprobability-Union{Tuple{A}, Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}, BioMarkovChain}} where A' href='#BioMarkovChains.dnaseqprobability-Union{Tuple{A}, Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}, BioMarkovChain}} where A'>#</a>&nbsp;<b><u>BioMarkovChains.dnaseqprobability</u></b> &mdash; <i>Method</i>.




```julia
dnaseqprobability(sequence::LongNucOrView{4}, model::BioMarkovChain)
```


Compute the probability of a given sequence using a transition probability matrix and the initial probabilities distributions of a `BioMarkovModel`.

$$P(X_1 = i_1, \ldots, X_T = i_T) = \pi_{i_1}^{T-1} \prod_{t=1}^{T-1} a_{i_t, i_{t+1}}$$

**Arguments**
- `sequence::LongNucOrView{4}`: The input sequence of nucleotides.
  
- `model::BioMarkovChain` is the actual data structure composed of a `tpm::Matrix{Float64}` the transition probability matrix and `initials=Vector{Float64}` the initial state probabilities.
  

**Returns**
- `probability::Float64`: The probability of the input sequence given the model.
  

**Example**

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



[source](https://github.com/camilogarciabotero/BioMarkovChains.jl/)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='BioMarkovChains.initials-Union{Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}}, Tuple{A}} where A' href='#BioMarkovChains.initials-Union{Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}}, Tuple{A}} where A'>#</a>&nbsp;<b><u>BioMarkovChains.initials</u></b> &mdash; <i>Method</i>.




```julia
initials(sequence::SeqOrView{A}) where A
```


Calculate the estimated initial probabilities for a Markov chain based on a given sequence.

This function takes a sequence of states and calculates the estimated initial probabilities of each state in the sequence for a Markov chain. The initial probabilities are estimated by counting the occurrences of each state at the beginning of the sequence and normalizing the counts to sum up to 1.

$$\begin{align}
\pi{i} &= P(X_{i} = i),  i \in T  \\
\sum_{i=1}^{N} \pi_{i} &= 1
\end{align}$$

Now using the dinucleotides counts estimating the initials would follow:

$$\hat{\pi_{i}} = c_{i} \sum_{k} c_{k}$$

**Arguments**
- `sequence::SeqOrView{A}`: The sequence of states representing the Markov chain.
  

**Returns**

An `Vector{Flot64}` of estimated initial probabilities for each state in the sequence.


[source](https://github.com/camilogarciabotero/BioMarkovChains.jl/)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='BioMarkovChains.log_odds_ratio_matrix-Tuple{BioMarkovChain, BioMarkovChain}' href='#BioMarkovChains.log_odds_ratio_matrix-Tuple{BioMarkovChain, BioMarkovChain}'>#</a>&nbsp;<b><u>BioMarkovChains.log_odds_ratio_matrix</u></b> &mdash; <i>Method</i>.




```julia
log_odds_ratio_matrix(model1::BioMarkovChain, model2::BioMarkovChain)
```


Calculates the log-odds ratio between the transition probability matrices of two BioMarkovChain models.

$$\beta = \log \frac{P(x|\mathscr{m}_{1})}{P(x|\mathscr{m}_{2})}$$

Where $\mathscr{m}_{1}$ and $\mathscr{m}_{2}$ are the two models transition probability matrices.

**Arguments**
- `model1::BioMarkovChain`: The first BioMarkovChain model.
  
- `model2::BioMarkovChain`: The second BioMarkovChain model.
  


[source](https://github.com/camilogarciabotero/BioMarkovChains.jl/)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='BioMarkovChains.log_odds_ratio_matrix-Union{Tuple{A}, Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}, BioMarkovChain}} where A' href='#BioMarkovChains.log_odds_ratio_matrix-Union{Tuple{A}, Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}, BioMarkovChain}} where A'>#</a>&nbsp;<b><u>BioMarkovChains.log_odds_ratio_matrix</u></b> &mdash; <i>Method</i>.




```julia
log_odds_ratio_matrix(sequence::NucleicSeqOrView{A}, model::BioMarkovChain) where A
```


Calculates the log-odds ratio between the transition probability matrix of a given DNA sequence and a reference model.

**Arguments**
- `sequence::NucleicSeqOrView{A}`: A DNA, RNA sequence or view with a length of 4 nucleotides.
  
- `model::BioMarkovChain`: A reference BioMarkovChain model.
  

**Examples**

```julia
sequence = LongNucOrView{4}("ACGT")
model = BioMarkovChain(sequence)  # Provide appropriate initialization for BioMarkovChain
result = log_odds_ratio_matrix(sequence, model)
```



[source](https://github.com/camilogarciabotero/BioMarkovChains.jl/)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='BioMarkovChains.log_odds_ratio_score-Union{Tuple{A}, Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}, BioMarkovChain}} where A' href='#BioMarkovChains.log_odds_ratio_score-Union{Tuple{A}, Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}, BioMarkovChain}} where A'>#</a>&nbsp;<b><u>BioMarkovChains.log_odds_ratio_score</u></b> &mdash; <i>Method</i>.




```julia
log_odds_ratio_score(sequence::SeqOrView{A}, model::BioMarkovChain; b::Number = ℯ)
```


Compute the log odds ratio score between a given sequence and a BioMarkovChain model.

$$S(x) = \sum_{i=1}^{L} \beta_{x_{i}x} = \sum_{i=1} \log \frac{a^{\mathscr{m}_{1}}_{i-1} x_i}{a^{\mathscr{m}_{2}}_{i-1} x_i}$$

**Arguments**
- `sequence::SeqOrView{A}`: A sequence of elements of type `A`.
  
- `model::BioMarkovChain`: A BioMarkovChain model.
  
- `b::Number = ℯ`: The base of the logarithm used to compute the log odds ratio.
  

**Returns**

The log odds ratio score between the sequence and the model.

**Example**


[source](https://github.com/camilogarciabotero/BioMarkovChains.jl/)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='BioMarkovChains.perronfrobenius-Union{Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}}, Tuple{A}} where A' href='#BioMarkovChains.perronfrobenius-Union{Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}}, Tuple{A}} where A'>#</a>&nbsp;<b><u>BioMarkovChains.perronfrobenius</u></b> &mdash; <i>Method</i>.




```julia
perronfrobenius(sequence::SeqOrView{A}, n::Int64=1) where A
```


Compute the Perron-Frobenius matrix, a column-stochastic version of the transition probability matrix (TPM), for a given nucleotide sequence.

The Perron-Frobenius matrix captures the asymptotic probabilities of transitioning between nucleotides in the sequence over a specified number of steps `n`. It provides insight into the long-term behavior of a Markov chain or a dynamical system associated with the sequence.

**Arguments**
- `sequence::SeqOrView{A}`: A nucleotide sequence represented as a `NucleicSeqOrView{A}` object.
  
- `n::Int64=1`: The number of steps to consider for the transition probability matrix. Default is 1.
  

**Returns**

A copy of the Perron-Frobenius matrix. Each column of this matrix corresponds to the probabilities of transitioning from the current nucleotide state to all possible nucleotide states after `n` steps.

**Example**

```julia
sequence = LongSequence{DNAAlphabet{4}}("ACGTCGTCCACTACGACATCAGC")  # Replace with an actual nucleotide sequence
n = 2
pf = perronfrobenius(sequence, n)
```



[source](https://github.com/camilogarciabotero/BioMarkovChains.jl/)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='BioMarkovChains.transition_count_matrix-Union{Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}}, Tuple{A}} where A' href='#BioMarkovChains.transition_count_matrix-Union{Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}}, Tuple{A}} where A'>#</a>&nbsp;<b><u>BioMarkovChains.transition_count_matrix</u></b> &mdash; <i>Method</i>.




```julia
transition_count_matrix(sequence::LongSequence{DNAAlphabet{4}})
```


Compute the transition count matrix (TCM) of a given DNA sequence.

**Arguments**
- `sequence::LongSequence{DNAAlphabet{4}}`: a `LongSequence{DNAAlphabet{4}}` object representing the DNA sequence.
  

**Returns**

A `Matrix` object representing the transition count matrix of the sequence.

**Example**

```
seq = LongDNA{4}("AGCTAGCTAGCT")

tcm = transition_count_matrix(seq)

4×4 Matrix{Int64}:
 0  0  3  0
 0  0  0  3
 0  3  0  0
 2  0  0  0
```



[source](https://github.com/camilogarciabotero/BioMarkovChains.jl/)

</div>
<br>
<div style='border-width:1px; border-style:solid; border-color:black; padding: 1em; border-radius: 25px;'>
<a id='BioMarkovChains.transition_probability_matrix-Union{Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}}, Tuple{A}, Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}, Int64}} where A' href='#BioMarkovChains.transition_probability_matrix-Union{Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}}, Tuple{A}, Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}, Int64}} where A'>#</a>&nbsp;<b><u>BioMarkovChains.transition_probability_matrix</u></b> &mdash; <i>Method</i>.




```julia
transition_probability_matrix(sequence::LongSequence{DNAAlphabet{4}}, n::Int64=1)
```


Compute the transition probability matrix (TPM) of a given DNA sequence. Formally it construct $\hat{\mathscr{M}}$ where: 

$$\mathscr{m}_{ij} = P(X_t = j \mid X_{t-1} = i) = \frac{{P(X_{t-1} = i, X_t = j)}}{{P(X_{t-1} = i)}}$$

The transition matrices of DNA and Amino-Acids are arranged sorted and in row-wise matrices:

First the DNA matrix:

$$\mathscr{M}_{DNA} = \begin{bmatrix}
_{AA} & _{AC} & _{AG} & _{AT} \\
_{CA} & _{CC} & _{CG} & _{CT} \\
_{GA} & _{GC} & _{GG} & _{GT} \\
_{TA} & _{TC} & _{TG} & _{TT} \\
\end{bmatrix}$$

And then, the Aminoacids:

$$\mathscr{M}_{AA} = \begin{bmatrix}
_{AA} & _{AC} & _{AD} & \dots & _{AW} \\
_{CA} & _{CC} & _{CD} & \dots & _{CW} \\
_{DA} & _{DC} & _{DD} & \dots & _{DW} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
_{WA} & _{WC} & _{WD} & \dots & _{WW} \\
\end{bmatrix}$$

**Arguments**
- `sequence::LongNucOrView{4}`: a `LongNucOrView{4}` object representing the DNA sequence.
  
- `n::Int64=1`: The order of the Markov model. That is the $\hat{M}^{n}$
  

**Keywords**
- `extended_alphabet::Bool=false`: If true will pass the extended alphabet of DNA to search
  

**Returns**

A `Matrix` object representing the transition probability matrix of the sequence.

**Example**

```
seq = dna"AGCTAGCTAGCT"

tpm = transition_probability_matrix(seq)

4×4 Matrix{Float64}:
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0
 0.0  1.0  0.0  0.0
 1.0  0.0  0.0  0.0
```



[source](https://github.com/camilogarciabotero/BioMarkovChains.jl/)

</div>
<br>
