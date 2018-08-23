using TimeseriesCausality
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

# Function aliases
anishchenko1 = TimeseriesCausality.Systems.anishchenko1
henon2 = TimeseriesCausality.Systems.henon2
henon4 = TimeseriesCausality.Systems.henon4
logistic2 = TimeseriesCausality.Systems.logistic2
logistic3 = TimeseriesCausality.Systems.logistic3
var1 = TimeseriesCausality.Systems.var1
var1coupled = TimeseriesCausality.Systems.var1coupled

# Check that initialisation happens correctly.
@test isa(anishchenko1(), DiscreteDynamicalSystem)
@test isa(henon2(), DiscreteDynamicalSystem)
@test isa(henon4(), DiscreteDynamicalSystem)
@test isa(logistic2(), DiscreteDynamicalSystem)
@test isa(logistic3(), DiscreteDynamicalSystem)
@test isa(var1(), DiscreteDynamicalSystem)
@test isa(var1coupled(), DiscreteDynamicalSystem)


# Initialise all the systems and generate trajectories
@test isa(trajectory(anishchenko1(), 10), Dataset)
@test isa(trajectory(henon2(), 10), Dataset)
@test isa(trajectory(henon4(), 10), Dataset)
@test isa(trajectory(logistic2(), 10), Dataset)
@test isa(trajectory(logistic3(), 10), Dataset)
@test isa(trajectory(var1(), 10), Dataset)
@test isa(trajectory(var1coupled(), 10), Dataset)


#################
# Continuous maps
#################

# Function aliases
chuacircuits_driven = TimeseriesCausality.Systems.chuacircuits_driven
chuacircuit_nscroll_sine = TimeseriesCausality.Systems.chuacircuit_nscroll_sine
rossler_rossler = TimeseriesCausality.Systems.rossler_rossler
rossler_lorenz = TimeseriesCausality.Systems.rossler_lorenz

# Check that initialisation happens correctly.
@test isa(chuacircuit_nscroll_sine(), ContinuousDynamicalSystem)
@test isa(chuacircuits_driven(), ContinuousDynamicalSystem)
@test isa(rossler_rossler(), ContinuousDynamicalSystem)
@test isa(rossler_lorenz(), ContinuousDynamicalSystem)

# Initialise all the systems and generate trajectories
@test isa(trajectory(chuacircuit_nscroll_sine(), 10), Dataset)
@test isa(trajectory(chuacircuits_driven(), 10), Dataset)
@test isa(trajectory(rossler_rossler(), 10), Dataset)
@test isa(trajectory(rossler_lorenz(), 10), Dataset)
