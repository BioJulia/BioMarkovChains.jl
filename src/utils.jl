"""
    transitions(sequence::LongSequence)

Compute the transition counts of each pair in a given biological sequence sequence.

# Arguments
- `sequence::LongSequence{DNAAlphabet{4}}`: a `LongSequence{DNAAlphabet{4}}` object representing the DNA sequence.

# Returns
A dictionary with keys being `Dict{Tuple{DNA, DNA}, Int64}` objects representing
the dinucleotides, and values being the number of occurrences of each dinucleotide
in the sequence.

# Example
```
seq = dna"AGCTAGCTAGCT"

dinucleotides(seq)
Dict{Tuple{DNA, DNA}, Int64} with 4 entries:
  (DNA_C, DNA_T) => 3
  (DNA_G, DNA_C) => 3
  (DNA_T, DNA_A) => 2
  (DNA_A, DNA_G) => 3
```
"""
function transitions(sequence::SeqOrView{A}) where A
    b = @view sequence[2:end]
    return countmap(zip(sequence, b))
end

"""
    hasprematurestop(sequence::LongNucOrView{4})::Bool

Determine whether the `sequence` of type `LongSequence{DNAAlphabet{4}}` contains a premature stop codon.

Returns a boolean indicating whether the `sequence` has more than one stop codon.
"""
function hasprematurestop(sequence::LongNucOrView{4})::Bool
    
    stopcodons = [LongDNA{4}("TAA"), LongDNA{4}("TAG"), LongDNA{4}("TGA")]  # Create a set of stop codons
    
    length(sequence) % 3 == 0 || error("The sequence is not divisible by 3")
    
    occursin(biore"T(AG|AA|GA)"dna, sequence[end-2:end]) || error("There is no stop codon at the end of the sequence")

    @inbounds for i in 1:3:length(sequence) - 4
        codon = sequence[i:i+2]
        if codon in stopcodons
            return true
        end
    end

    return false
end

function _int_to_dna(index::Int64; extended_alphabet::Bool = false)
    A = extended_alphabet ? collect(alphabet(DNA)) : [DNA_A, DNA_C, DNA_G, DNA_T]
    return LongSequence{DNAAlphabet{4}}([A[index]])
end

function _dna_to_int(nucleotide::DNA; extended_alphabet::Bool = false)
    A = extended_alphabet ? collect(alphabet(DNA)) : [DNA_A, DNA_C, DNA_G, DNA_T]
    return findfirst(nucleotide, LongSequence{DNAAlphabet{4}}(A))
end

# function randbmc(A::DataType, n::Int64=1)

#     if A == DNA || A == RNA
#         tpm = rand(4,4)
#         inits = rand(4)
#     elseif A == AminoAcid
#         tpm = rand(20,20)
#         inits = rand(20)
#     else
#         error("Alphabet must be of the DNA, RNA or AminoAcid DataType")
#     end

#     BMC(tpm, inits, n)
# end

function randbmc(A::DataType, n::Int64=1)
    if A == DNA || A == RNA
        tpm = rand(4, 4)
        row_sums = sum(tpm, dims=2)
        inits = rand(4)
        init_sum = sum(inits)
        @views tpm ./= row_sums
        @views inits ./= init_sum
    elseif A == AminoAcid
        tpm = rand(20, 20)
        row_sums = sum(tpm, dims=2)
        inits = rand(20)
        init_sum = sum(inits)
        @views tpm ./= row_sums
        @views inits ./= init_sum
    else
        error("Alphabet must be of the DNA, RNA or AminoAcid DataType")
    end

    BMC(tpm, inits, n)
end

