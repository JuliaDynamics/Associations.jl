using Test
using CausalityTools
using DelayEmbeddings, DynamicalSystemsBase

include("methods/test_smeasure.jl")
include("methods/test_joint_distance_distribution.jl")
include("methods/test_predictive_asymmetry.jl")
include("methods/test_crossmapping.jl")
include("methods/test_convergentcrossmapping.jl")

include("systems/discrete/test_discrete_systems.jl")
include("systems/continuous/test_continuous_systems.jl")

#include("integrations/test_uncertaindata_integration.jl")
