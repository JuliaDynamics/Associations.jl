using CausalityTools
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

using DynamicalSystems
import DynamicalSystemsBase:
    DynamicalSystem,
    DiscreteDynamicalSystem,
    ContinuousDynamicalSystem,
    Dataset

###############
# Discrete maps
###############

# Check that initialisation happens correctly.
@test isa(ar1(), DiscreteDynamicalSystem)
@test isa(anishchenko1(), DiscreteDynamicalSystem)
@test isa(henon2(), DiscreteDynamicalSystem)
@test isa(henon4(), DiscreteDynamicalSystem)
@test isa(logistic2(), DiscreteDynamicalSystem)
@test isa(logistic3(), DiscreteDynamicalSystem)
@test isa(logistic4(), DiscreteDynamicalSystem)
@test isa(var1(), DiscreteDynamicalSystem)
@test isa(var1coupled(), DiscreteDynamicalSystem)
@test isa(verdes(), DiscreteDynamicalSystem)

# Discrete maps with lags directly return a Dataset
@test isa(henon_triple(n = 10), Dataset)

# Initialise all the systems and generate trajectories
@test isa(trajectory(ar1(), 10), Dataset)
@test isa(trajectory(anishchenko1(), 10), Dataset)
@test isa(trajectory(henon2(), 10), Dataset)
@test isa(trajectory(henon4(), 10), Dataset)
@test isa(trajectory(logistic2(), 10), Dataset)
@test isa(trajectory(logistic3(), 10), Dataset)
@test isa(trajectory(logistic4(), 10), Dataset)
@test isa(trajectory(var1(), 10), Dataset)
@test isa(trajectory(var1coupled(), 10), Dataset)
@test isa(trajectory(verdes(), 10), Dataset)


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
