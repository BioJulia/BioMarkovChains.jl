import Base: show

function Base.show(io::IO, tcm::TransitionCountMatrix)
    nucleotides = sort(collect(keys(tcm.order)))

    # Print type
    println(io, "TCM{$(typeof(tcm.order)), $(typeof(tcm.counts)):")

    # Print header
    max_digits = maximum([length(string(maximum(tcm.counts[:, i]))) for i in 1:size(tcm.counts, 2)])
    header_str = "   " * join([rpad(n, max_digits+1) for n in nucleotides], "")
    println(io, header_str)

    # Print rows
    for (i, nucleotide1) in enumerate(nucleotides)
        row = tcm.counts[i, :]
        row_str = join([rpad(string(x), max_digits+1) for x in row], "")
        println(io, "$nucleotide1  $row_str")
    end
end

function Base.show(io::IO, tpm::TransitionProbabilityMatrix)
    nucleotides = sort(collect(keys(tpm.order)))

    # Print type
    println(io, "TPM{$(typeof(tpm.order)), $(typeof(tpm.probabilities)):")

    # Print header
    max_digits = maximum([length(string(maximum(round.(tpm.probabilities[:, i], digits = 3)))) for i in 1:size(tpm.probabilities, 2)])
    header_str = "   " * join([rpad(n, max_digits+1) for n in nucleotides], "")
    println(io, header_str)

    # Print rows
    for (i, nucleotide1) in enumerate(nucleotides)
        row = round.(tpm.probabilities[i, :], digits = 3)
        row_str = join([rpad(string(x), max_digits+1) for x in row], "")
        println(io, "$nucleotide1  $row_str")
    end
end

function Base.show(io::IO, model::TransitionModel)
    # Print the type name
    println(io, "TransitionModel:")

    # Print the size of the transition probability matrix
    tpm_size = "  - Transition Probability Matrix (Size: $(size(model.tpm, 1)) × $(size(model.tpm, 2))):"
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
    initials_size = "  - Initials (Size: $(size(model.initials, 1)) × $(size(model.initials, 2))):"
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
    order = "  - order: $(model.n)"
    println(io, order)
end

Base.length(::TransitionModel) = 1