
####################
# Continuous systems
####################

# Check that initialisation happens correctly.
@test isa(chuacircuit_nscroll_sine(), ContinuousDynamicalSystem)
#@test isa(chuacircuits_driven(), ContinuousDynamicalSystem)
@test isa(hindmarsh_rose(), ContinuousDynamicalSystem)
@test isa(lorenz_triple(), ContinuousDynamicalSystem)
@test isa(mediated_link(), ContinuousDynamicalSystem)
@test isa(rossler_rossler(), ContinuousDynamicalSystem)
@test isa(rossler_lorenz(), ContinuousDynamicalSystem)

# Initialise all the systems and generate trajectories
l = 0.1
@test isa(trajectory(chuacircuit_nscroll_sine(), l), Dataset)
#@test isa(trajectory(chuacircuits_driven(), l), Dataset)
@test isa(trajectory(hindmarsh_rose(), l), Dataset)
@test isa(trajectory(lorenz_triple(), l), Dataset)
@test isa(trajectory(mediated_link(), l), Dataset)
@test isa(trajectory(rossler_rossler(), l), Dataset)
@test isa(trajectory(rossler_lorenz(), l), Dataset)
