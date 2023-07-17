import Base: show

function Base.show(io::IO, model::TransitionModel)
    # Print the type name
    println(io, "TransitionModel:")

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
    initials_size = "  - Initial Probabilities -> Vector{Float64}($(size(model.initials, 1)) × $(size(model.initials, 2))):"
    println(io, initials_size)
    for row in 1:size(model.initials, 1)
        print(io, "    ")
        for col in 1:size(model.initials, 2)
            print(io, round(model.initials[row, col], digits=3))
            print(io, "\t")
        end
        println(io)
    end

    # Print the value of 'n'
    order = "  - Markov Chain Order:$(model.n)"
    println(io, order)
end

Base.length(::TransitionModel) = 1