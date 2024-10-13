@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    seq = LongDNA{4}("ACTACTACTACTACTA")
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        transition_count_matrix(seq)
        transition_probability_matrix(seq)

    end
end