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

##################
# Function aliases
##################
logistic3 = TimeseriesCausality.Systems.logistic3


###############
# Discrete maps
###############
@test isa(logistic3(), DiscreteDynamicalSystem)


############################################################
# Initialise all the systems and generate trajectories
############################################################
@test isa(trajectory(logistic3(), 10), Dataset)
