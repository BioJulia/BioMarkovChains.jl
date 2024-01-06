import Base: show, length, eltype, size
# import StatsAPI: fit!

function Base.show(
    io::IO, 
    model::BioMarkovChain
)
    # # Print the type name
    # println(io, "BioMarkovChain:")

    # Determine the alphabet type
    # alphabet_type = length(model) == 20 ? "AminoAcids" : "DNA/RNA"
    alphabet_type = eltype(model)

    # Print the type name with inferred alphabet type
    println(io, "BioMarkovChain of $alphabet_type:")

    # Print the transition probability matrix
    println(io, "  - Transition Probability Matrix -> Matrix{Float64}($(size(model.tpm, 1)) × $(size(model.tpm, 2))):")
    for row in 1:size(model.tpm, 1)
        print(io, "  ")
        max_cols = size(model.tpm, 2)
        cols_to_show = min(6, max_cols)
        
        for col in 1:cols_to_show
            print(io, " ", rpad(round(model.tpm[row, col], digits=4), 7))
        end
        
        if max_cols > 6
            print(io, "   ", rpad("…", 5))
            
            for col in (max_cols - 6):max_cols
                print(io, " ", rpad(round(model.tpm[row, col], digits=4), 7))
            end
        end
        
        println(io)
    end

    # Print the initials matrix
    println(io, "  - Initial Probabilities -> Vector{Float64}($(size(model.inits, 1)) × $(size(model.inits, 2))):")
    for row in 1:size(model.inits, 1)
        println(io, "   ", round(model.inits[row], digits=4))
    end

    # Print the value of 'n'
    println(io, "  - Markov Chain Order -> Int64:")
    println(io, "   ", "$(model.n)")
end

@inline Base.length(bmc::BioMarkovChain) = length(bmc.inits)
@inline Base.size(bmc::BioMarkovChain) = size(bmc.tpm)
@inline Base.eltype(bmc::BioMarkovChain) = bmc.alphabet

"""
    fit!(bmc::BMC, inits:Vector{Float64}, tpm::Matrix{Float64})

Update `bmc` in-place based on information generated from a state sequence.
"""
# function StatsAPI.fit!(bmc::BMC, inits::Vector{Float64}, tpm::Matrix{Float64})
#     bmc.inits .= inits
#     sum_to_one!(bmc.inits)
#     bmc.tpm .= tpm
#     foreach(sum_to_one!, eachrow(bmc.tpm))
#     return nothing
# end