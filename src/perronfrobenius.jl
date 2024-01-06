### MarkdovChainHammer.jl ###

"""
    perronfrobenius(sequence::SeqOrView{A}, n::Int64=1) where A

Compute the Perron-Frobenius matrix, a column-stochastic version of the transition probability matrix (TPM), for a given nucleotide sequence.

The Perron-Frobenius matrix captures the asymptotic probabilities of transitioning between nucleotides in the sequence over a specified number of steps `n`. It provides insight into the long-term behavior of a Markov chain or a dynamical system associated with the sequence.

# Arguments
- `sequence::SeqOrView{A}`: A nucleotide sequence represented as a `NucleicSeqOrView{A}` object.
- `n::Int64=1`: The number of steps to consider for the transition probability matrix. Default is 1.

# Returns
A copy of the Perron-Frobenius matrix. Each column of this matrix corresponds to the probabilities of transitioning from the current nucleotide state to all possible nucleotide states after `n` steps.

# Example
```julia
sequence = LongSequence{DNAAlphabet{4}}("ACGTCGTCCACTACGACATCAGC")  # Replace with an actual nucleotide sequence
n = 2
pf = perronfrobenius(sequence, n)
```
"""
function perronfrobenius(
    sequence::SeqOrView{A};
    n::Int64=1
) where {A}
    tpm = transition_probability_matrix(sequence, n)
    return copy(tpm')
end

function perronfrobenius(
    bmc::BioMarkovChain
)
   return copy(bmc.tpm')
end

function perronfrobenius(
    tpm::Matrix{Float64}
)
    return copy(tpm')
end

"""
    generatedna(pf::Matrix{Float64}, steps::Int64; extended_alphabet::Bool = false)

Generate a new DNA sequence using the Perron-Frobenius matrix-based approach.

This function utilizes the Perron-Frobenius matrix `pf` to generate a DNA sequence over a specified number of steps. Each step involves transitioning from the current DNA state to a new state based on the probabilities defined in the matrix `pf`.

# Arguments
- `pf::Matrix{Float64}`: The Perron-Frobenius matrix, a square matrix of transition probabilities. Each row corresponds to the source state, and each column corresponds to the probabilities of transitioning to different states.
- `steps::Int64`: The number of steps to generate the new DNA sequence.
- `extended_alphabet::Bool = false`: (Optional) If `true`, the function will use an extended DNA alphabet with additional characters beyond the standard 'A', 'C', 'G', and 'T'. Default is `false`.

# Returns
A `LongDNA{4}` object representing the newly generated DNA sequence.

# Example
```julia
using GeneFinder, BioSequences

seq = randdnaseq(10^6) # create a random DNA sequence of 1Mb
orfdna = getorfdna(seq, min_len=90)[1] # find a random ORF

pf = perronfrobenius(orfdna)

    4×4 Matrix{Float64}:
    0.333333  0.296296  0.238095  0.292683
    0.245614  0.314815  0.238095  0.317073
    0.210526  0.222222  0.238095  0.195122
    0.210526  0.166667  0.285714  0.195122

newdna = generatedna(pf, 100)
```
"""
# function generatedna(pf::Matrix{Float64}, steps::Int64; extended_alphabet::Bool = false)
#     newseq = LongDNA{4}()
#     # pf = transpose(tpm) # The Perron-Frobenius matrix
#     trajectory = generate(pf, steps)
#     for i in trajectory
#         newseq = append!(newseq, _int_to_dna(i; extended_alphabet))
#     end
#     return newseq
# end