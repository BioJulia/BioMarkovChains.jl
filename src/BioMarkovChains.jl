module BioMarkovChains

using BioSequences:
    LongSequence,
    LongSubSeq,
    LongNuc,
    LongDNA,
    LongRNA,
    ACGT,
    ACGU,
    NucleicAcidAlphabet,
    DNA,
    DNA_A, DNA_C, DNA_G, DNA_T, DNA_M,
    DNAAlphabet,
    RNAAlphabet,
    Alphabet,
    
    #RNA
    RNA,
    
    #AminoAcids
    AminoAcid,
    AminoAcidAlphabet,
    LongAA,
    AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V, AA_O, AA_U, AA_B, AA_J, AA_Z, AA_X, AA_Term, AA_Gap,

    # Other functions
    SeqOrView, NucleicSeqOrView

    #tests and precompilation

using PrecompileTools: @setup_workload, @compile_workload
using VectorizedKmers: count_kmers

include("types.jl")
include("utils.jl")
include("transitions.jl")
include("models.jl")
include("perronfrobenius.jl")
include("extended.jl")

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    seq = LongDNA{4}("ACTACTACTACTACTA")
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        transition_count_matrix(seq)
        transition_probability_matrix(seq)

    end
end

end # BioMarkovChains
