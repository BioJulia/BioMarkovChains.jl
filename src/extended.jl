import Base: show, length, size
import BioSequences: Alphabet
# import StatsAPI: fit!

function Base.show(
    io::IO, 
    model::BioMarkovChain{A}
) where {A}

    # Determine the alphabet type
    alphtype = eltype(model)

    # Print the type name with inferred alphabet type
    println(io, "BioMarkovChain of $alphtype alphabet and order $(model.n):")

    # Print the transition probability matrix
    println(io, "  - Transition Probability Matrix -> Matrix{Float64}($(size(model.tpm, 1)) × $(size(model.tpm, 2))):")
    for row in 1:size(model.tpm, 1)
        print(io, "  ")
        max_cols = size(model.tpm, 2)
        cols_to_show = min(6, max_cols)
        
        for col in 1:cols_to_show
            print(io, " ", rpad(round(model.tpm[row, col], digits=3), 6))
        end
        
        if max_cols > 6
            print(io, "   ", rpad("…", 4))
            
            for col in (max_cols - 6):max_cols
                print(io, " ", rpad(round(model.tpm[row, col], digits=3), 6))
            end
        end
        
        println(io)
    end

    # Print the initials matrix
    println(io, "  - Initial Probabilities -> Vector{Float64}($(size(model.inits, 1)) × $(size(model.inits, 2))):")
    print(io, "  ")
    
    max_initials = length(model.inits)
    initials_to_show = min(6, max_initials)
    
    for row in 1:initials_to_show
        print(io, " ", rpad(round(model.inits[row], digits=3), 6))
    end
    
    if max_initials > 6
        print(io, "   ", rpad("…", 4))
        
        for row in (max_initials - 6):max_initials
            print(io, " ", rpad(round(model.inits[row], digits=3), 6))
        end
    end
    
end

@inline Base.length(bmc::BioMarkovChain{A}) where {A} = length(bmc.inits)
@inline Base.size(bmc::BioMarkovChain{A}) where {A} = size(bmc.tpm)
@inline Base.eltype(bmc::BioMarkovChain{A}) where {A} = eltype(A)

## Overload operators

Base.:(==)(a::BioMarkovChain{A}, b::BioMarkovChain{A}) where {A} = eltype(a) == eltype(b) && a.tpm == b.tpm && a.inits == b.inits
Base.:(^)(a::BioMarkovChain{A}, n::Int) where {A} = BioMarkovChain{A}(a.tpm^n, a.inits, n)

## BioSequences

Alphabet(bmc::BioMarkovChain{A}) where {A} = A()
# Alphabet(::Type{BioMarkovChain{A}}) where {A} = A

## Random

# mutable struct BioMarkovChainDist{T}
#     bmcmean::T
# end

# function Random.rand(rng::AbstractRNG, dist::BioMarkovChainDist)
#     quantity = dist.bmcmean + randn(rng)
#     return Stuff(quantity)
# end

# using Random: Random, randn
# import Random: Random, randn, AbstractRNG

# export Stuff, StuffDist
# struct Stuff{T}
#     quantity::T
# end

# mutable struct StuffDist{T}
#     quantity_mean::T
# end

# function Random.rand(rng::AbstractRNG, dist::StuffDist)
#     quantity = dist.quantity_mean + randn(rng)
#     return Stuff(quantity)
# end
