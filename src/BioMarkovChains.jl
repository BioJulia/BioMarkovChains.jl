module BioMarkovChains

using BioSequences:
    ACGT,
    DNA,
    DNA_A,
    DNA_C,
    DNA_G,
    DNA_T,
    DNAAlphabet,
    NucleicAcidAlphabet,
    LongSequence,
    LongSubSeq,
    LongDNA,
    alphabet,
    @biore_str

using MarkovChainHammer.Trajectory: generate
using PrecompileTools: @setup_workload, @compile_workload
using TestItems: @testitem

include("types.jl")
export TransitionCountMatrix, TransitionProbabilityMatrix, TransitionModel

include("utils.jl")
export dinucleotides, hasprematurestop

include("transitions.jl")
export transition_count_matrix, transition_probability_matrix, transition_model, initial_distribution, sequenceprobability

include("models.jl")
export ECOLICDS, ECOLINOCDS

include("perronfrobenius.jl")
export perronfrobenius, generatedna

include("extended.jl")

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    using BioSequences
    seq = randdnaseq(10^3)
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        dinucleotides(seq)
        sequenceprobability(seq, ECOLICDS)
        transition_model(seq)
        perronfrobenius(seq)
    end
end

end # BioMarkovChains
