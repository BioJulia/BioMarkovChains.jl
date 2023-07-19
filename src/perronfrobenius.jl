### MarkdovChainHammer.jl ###

function perronfrobenius(sequence::LongNucOrView{4}, n::Int64=1)
    tpm = transition_probability_matrix(sequence, n)
    pf = transpose(tpm)
    return copy(pf)
end

function perronfrobenius(bmc::BioMarkovChain)
   tpm = bmc.tpm
   pf = transpose(tpm)
   return copy(pf)
end

function perronfrobenius(tpm::Matrix{Float64})
    pf = transpose(tpm)
    return copy(pf)
end

function generatedna(pf::Matrix{Float64}, steps::Int64; extended_alphabet::Bool = false)
    newseq = LongDNA{4}()
    # pf = transpose(tpm) # The Perron-Frobenius matrix
    trajectory = generate(pf, steps)
    for i in trajectory
        newseq = append!(newseq, _int_to_dna(i; extended_alphabet))
    end
    return newseq
end