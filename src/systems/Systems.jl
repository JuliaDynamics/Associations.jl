module Systems

using Distributions
using DynamicalSystems
import DynamicalSystemsBase:
    DynamicalSystem,
    ContinuousDynamicalSystem,
    DiscreteDynamicalSystem,
    Dataset

import DifferentialEquations: @ode_def

include("discretemaps/logistic3.jl")


# We export both the equations of motion and function to generate
# the systems for every example. The equations of motion functions
# are always prepended with `eom_`.
export
###############
# Discrete maps
###############
eom_logistic3, logistic3


############################################################
# Initialise all the systems once, generating a trajectory.
############################################################
trajectory(logistic3(), 1)

end
