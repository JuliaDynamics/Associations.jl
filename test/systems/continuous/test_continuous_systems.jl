
####################
# Continuous systems
####################
@testset "Continuous example systems" begin 

    # Check that initialisation happens correctly.
    #@test isa(chuacircuits_driven(), ContinuousDynamicalSystem)
    @test isa(chuacircuit_nscroll_sine(), ContinuousDynamicalSystem)
    @test isa(hindmarsh_rose(), ContinuousDynamicalSystem)

    @test isa(lorenz_lorenz_bidir(), ContinuousDynamicalSystem)
    @test isa(lorenz_lorenz_lorenz_bidir_forced(), ContinuousDynamicalSystem)
    @test isa(lorenz_lorenz_lorenz_transitive(), ContinuousDynamicalSystem)
    @test isa(rossler_rossler_bidir(), ContinuousDynamicalSystem)
    @test isa(rossler_rossler_rossler_bidir_forced(), ContinuousDynamicalSystem)

    @test isa(mediated_link(), ContinuousDynamicalSystem)

    @test isa(rossler_lorenz(), ContinuousDynamicalSystem)

    # Initialise all the systems and generate trajectories
    l = 0.1
    @test isa(trajectory(chuacircuit_nscroll_sine(), l), Dataset)
    #@test isa(trajectory(chuacircuits_driven(), l), Dataset)
    @test isa(trajectory(hindmarsh_rose(), l), Dataset)

    @test isa(trajectory(lorenz_lorenz_bidir(), l), Dataset)
    @test isa(trajectory(lorenz_lorenz_lorenz_bidir_forced(), l), Dataset)
    @test isa(trajectory(lorenz_lorenz_lorenz_transitive(), l), Dataset)

    @test isa(trajectory(rossler_rossler_bidir(), l), Dataset)
    @test isa(trajectory(rossler_rossler_rossler_bidir_forced(), l), Dataset)

    @test isa(trajectory(mediated_link(), l), Dataset)

    @test isa(trajectory(rossler_lorenz(), l), Dataset)

end