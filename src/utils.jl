function _int_to_dna(index::Int64; extended_alphabet::Bool = false)
    A = extended_alphabet ? alphabet(DNA) : ACGT
    return LongSequence{DNAAlphabet{4}}([A[index]])
end

function _dna_to_int(nucleotide::DNA; extended_alphabet::Bool = false)
    A = extended_alphabet ? [DNA_Gap, DNA_A, DNA_C, DNA_M, DNA_G, DNA_R, DNA_S, DNA_V, DNA_T, DNA_W, DNA_Y, DNA_H, DNA_K, DNA_D, DNA_B, DNA_N] : [DNA_A, DNA_C, DNA_G, DNA_T]
    return searchsortedfirst(A, nucleotide) # findfirst(nucleotide, LongSequence{DNAAlphabet{4}}(A))
end

function randbmc(statespace::DataType, n::Int64=1)
    if statespace == DNA || statespace == RNA
        tpm = rand(4, 4)
        row_sums = sum(tpm, dims=2)
        inits = rand(4)
        init_sum = sum(inits)
        @views tpm ./= row_sums
        @views inits ./= init_sum
    elseif statespace == AminoAcid
        tpm = rand(20, 20)
        row_sums = sum(tpm, dims=2)
        inits = rand(20)
        init_sum = sum(inits)
        @views tpm ./= row_sums
        @views inits ./= init_sum
    else
        error("Alphabet must be of the DNA, RNA or AminoAcid DataType")
    end

    return BMC(statespace, tpm, inits, n)
end