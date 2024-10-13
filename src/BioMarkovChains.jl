module BioMarkovChains

using BioSequences:
    # LongSequence,
    LongDNA,
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
include("workload.jl")

end # BioMarkovChains
