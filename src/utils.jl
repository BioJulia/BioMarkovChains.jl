function _int_to_dna(index::Int64)
    # A = extended_alphabet ? alphabet(DNA) : ACGT
    modifier(value) = (value == 3) ? 4 : (value == 4) ? 8 : value
    return reinterpret.(DNA, Int8(modifier(index))) #LongSequence{DNAAlphabet{4}}([A[index]]) # reinterpret.(DNA, Int8(1)) == A, reinterpret.(Int8, DNA_A) = 1
end

function _dna_to_int(nucleotide::DNA)
    # A = extended_alphabet ? [DNA_Gap, DNA_A, DNA_C, DNA_M, DNA_G, DNA_R, DNA_S, DNA_V, DNA_T, DNA_W, DNA_Y, DNA_H, DNA_K, DNA_D, DNA_B, DNA_N] : [DNA_A, DNA_C, DNA_G, DNA_T]
    modifier(value) = (value == DNA_G) ? DNA_M : (value == DNA_T) ? DNA_G : value
    return reinterpret.(Int8, modifier(nucleotide))[1]  #searchsortedfirst(A, nucleotide) # findfirst(nucleotide, LongSequence{DNAAlphabet{4}}(A))
end

function randbmc(A::Alphabet, n::Int64=1)

    if !(A in (DNAAlphabet{4}(), RNAAlphabet{4}(), AminoAcidAlphabet()))
        throw(ArgumentError("Alphabet must be of the DNAAlphabet, RNAAlphabet, or AminoAcidAlphabet."))
    end

    nstates = (A == AminoAcidAlphabet) ? 20 : 4
    tpm = rand(nstates, nstates)
    
    # Normalize rows of the transition probability matrix
    rowsums = sum(tpm, dims=2)
    @views tpm ./= rowsums
    
    # Generate random initial probabilities and normalize them
    inits = rand(nstates)
    initsum = sum(inits)
    @views inits ./= initsum
    
    return BMC(A, tpm, inits, n)
end
