# From BMC to sequence probabilties

We can now calculate a transition probability matrix from a `LongDNA` sequence using the `transition_probability_matrix` and `initials` methods:

``` julia
using BioSequences

sequence = dna"CCTCCCGGACCCTGGGCTCGGGAC"

tpm = transition_probability_matrix(sequence)
initials = initials(sequence)

println(tpm)
println(initials)
```

    4×4 Matrix{Float64}:
     0.0   1.0       0.0       0.0
     0.0   0.5       0.2       0.3
     0.25  0.125     0.625     0.0
     0.0   0.666667  0.333333  0.0

    4-element Vector{Float64}:
     0.08695652173913043
     0.43478260869565216
     0.34782608695652173
     0.13043478260869565
```

