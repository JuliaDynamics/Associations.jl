
####################
# Continuous systems
####################

# Check that initialisation happens correctly.
@test isa(chuacircuit_nscroll_sine(), ContinuousDynamicalSystem)
@test isa(chuacircuits_driven(), ContinuousDynamicalSystem)
@test isa(hindmarsh_rose(), ContinuousDynamicalSystem)
@test isa(lorenz_triple(), ContinuousDynamicalSystem)
@test isa(mediated_link(), ContinuousDynamicalSystem)
@test isa(rossler_rossler(), ContinuousDynamicalSystem)
@test isa(rossler_lorenz(), ContinuousDynamicalSystem)

# Initialise all the systems and generate trajectories
@test isa(trajectory(chuacircuit_nscroll_sine(), 10), Dataset)
@test isa(trajectory(chuacircuits_driven(), 10), Dataset)
@test isa(trajectory(hindmarsh_rose(), 10), Dataset)
@test isa(trajectory(lorenz_triple(), 10), Dataset)
@test isa(trajectory(mediated_link(), 10), Dataset)
@test isa(trajectory(rossler_rossler(), 10), Dataset)
@test isa(trajectory(rossler_lorenz(), 10), Dataset)
