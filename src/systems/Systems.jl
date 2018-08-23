module Systems

using Distributions
using DynamicalSystems
import DynamicalSystemsBase:
    DynamicalSystem,
    ContinuousDynamicalSystem,
    DiscreteDynamicalSystem,
    Dataset

import DifferentialEquations: @ode_def

include("discretemaps/henon2.jl")
include("discretemaps/logistic2.jl")
include("discretemaps/logistic3.jl")
include("discretemaps/var1.jl")


# We export both the equations of motion and function to generate
# the systems for every example. The equations of motion functions
# are always prepended with `eom_`.
export
###############
# Discrete maps
###############
eom_henon2, henon2
eom_logistic2, logistic2
eom_logistic3, logistic3
eom_var1, var1

# Initialise all the systems once, generating a trajectory.
trajectory(henon2(), 10)
trajectory(logistic2(), 10)
trajectory(logistic3(), 10)
trajectory(var1(), 10)

####################
# Continuous systems
###################
eom_rossler_rossler

# Initialise all the systems once, generating a trajectory.
trajectory(eom_rossler_rossler(), 10)



end
