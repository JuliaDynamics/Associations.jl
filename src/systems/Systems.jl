module Systems

using Distributions
using DynamicalSystems
import DynamicalSystemsBase:
    DynamicalSystem,
    ContinuousDynamicalSystem,
    DiscreteDynamicalSystem,
    Dataset

import DifferentialEquations: @ode_def

include("discretemaps/anishchenko1.jl")
include("discretemaps/henon2.jl")
include("discretemaps/henon4.jl")
include("discretemaps/logistic2.jl")
include("discretemaps/logistic3.jl")
include("discretemaps/var1.jl")
include("discretemaps/var1coupled.jl")

include("continuous_systems/chuacircuits_driven.jl")
include("continuous_systems/chuacircuit_nscroll_sine.jl")
include("continuous_systems/rosslerrossler.jl")
include("continuous_systems/rosslerlorenz.jl")

# We export both the equations of motion and function to generate
# the systems for every example. The equations of motion functions
# are always prepended with `eom_`.
export
###############
# Discrete maps
###############
eom_anishchenko1, anishchenko1,
eom_henon2, henon2,
eom_henon4, henon4,
eom_logistic2, logistic2,
eom_logistic3, logistic3,
eom_var1, var1,
eom_var1coupled, var1coupled


export
####################
# Continuous systems
###################
eom_chuacircuit_nscroll_sine, chuacircuit_nscroll_sine,
eom_chuacircuits_driven, chuacircuits_driven,
eom_rossler_rossler, rossler_rossler,
eom_rossler_lorenz, rossler_lorenz

# Initialise all the systems once, generating a trajectory.
trajectory(anishchenko1(), 10)
trajectory(henon2(), 10)
trajectory(henon4(), 10)
trajectory(logistic2(), 10)
trajectory(logistic3(), 10)
trajectory(var1(), 10)
trajectory(var1coupled(), 10)

trajectory(rossler_rossler(), 10)
trajectory(rossler_lorenz(), 10)
trajectory(chuacircuits_driven(), 10)
trajectory(chuacircuit_nscroll_sine(), 10)


end
