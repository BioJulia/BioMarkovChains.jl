import Base: show, length

function Base.show(io::IO, model::BioMarkovChain)
    # Print the type name
    println(io, "BioMarkovChain:")

    # Print the size of the transition probability matrix
    tpm_size = "  - Transition Probability Matrix -> Matrix{Float64}($(size(model.tpm, 1)) × $(size(model.tpm, 2))):"
    println(io, tpm_size)
    for row in 1:size(model.tpm, 1)
        print(io, "    ")
        for col in 1:size(model.tpm, 2)
            print(io, round(model.tpm[row, col], digits=3))
            print(io, "\t")
        end
        println(io)
    end

    # Print the size of the initials matrix
    initials_size = "  - Initial Probabilities -> Vector{Float64}($(size(model.inits, 1)) × $(size(model.inits, 2))):"
    println(io, initials_size)
    for row in 1:size(model.inits, 1)
        print(io, "    ")
        for col in 1:size(model.inits, 2)
            print(io, round(model.inits[row, col], digits=3))
            print(io, "\t")
        end
        println(io)
    end

    # Print the value of 'n'
    order = "  - Markov Chain Order:$(model.n)"
    println(io, order)
end

Base.length(bmc::BioMarkovChain) = length(bmc.inits)