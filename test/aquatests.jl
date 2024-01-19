@testset "Aqua" begin

    Aqua.test_ambiguities(BioMarkovChains)
    Aqua.test_persistent_tasks(BioMarkovChains)
    Aqua.test_piracies(BioMarkovChains)
    Aqua.test_stale_deps(BioMarkovChains)
    Aqua.test_unbound_args(BioMarkovChains)
    Aqua.test_undefined_exports(BioMarkovChains)

end