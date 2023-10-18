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
    DNAAlphabet,
    DNA_Gap, DNA_A, DNA_C, DNA_M, DNA_G, DNA_R, DNA_S, DNA_V, DNA_T, DNA_W, DNA_Y, DNA_H, DNA_K, DNA_D, DNA_B, DNA_N,
    
    #RNA
    RNA,
    RNA_Gap, RNA_A, RNA_C, RNA_M, RNA_G, RNA_R, RNA_S, RNA_V, RNA_U, RNA_W, RNA_Y, RNA_H, RNA_K, RNA_D, RNA_B, RNA_N,
    
    #AminoAcids
    AminoAcid,
    AminoAcidAlphabet,
    LongAA,
    AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V, AA_O, AA_U, AA_B, AA_J, AA_Z, AA_X, AA_Term, AA_Gap,

    # Other functions
    SeqOrView,
    alphabet,
    @biore_str
    

using MarkovChainHammer.Trajectory: generate
using PrecompileTools: @setup_workload, @compile_workload
using TestItems: @testitem
using StatsAPI: StatsAPI, fit, fit!
using StatsBase: countmap
using VectorizedKmers: count_kmers

include("types.jl")
export BioMarkovChain, BMC

include("utils.jl")
export transitions, hasprematurestop, randbmc

include("transitions.jl")
export initials, transition_count_matrix, transition_probability_matrix, dnaseqprobability, logoddsratio, logoddsratioscore

include("models.jl")
export ECOLICDS, ECOLINOCDS

include("perronfrobenius.jl")
export perronfrobenius, generatedna

include("extended.jl")

include("constants.jl")

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    using BioSequences
    seq = randdnaseq(10^3)
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        transitions(seq)
        transition_count_matrix(seq)
        transition_probability_matrix(seq)
        dnaseqprobability(seq, ECOLICDS)
        BioMarkovChain(seq)
        perronfrobenius(seq)
    end
end

end # BioMarkovChains
