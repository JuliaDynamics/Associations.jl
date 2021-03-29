@reexport module ExampleSystems

using Distributions
using StatsBase
using DynamicalSystems

include("discretemaps/ar1.jl")
include("discretemaps/ar1_bidir.jl")
include("discretemaps/anishchenko1.jl")
include("discretemaps/henon2.jl")
include("discretemaps/henontriple.jl")
include("discretemaps/ikeda.jl")
include("discretemaps/linearmap.jl")
include("discretemaps/nonlinear3D_linear_and_nonlinear_coupling.jl")
include("discretemaps/nontrivial_pegiun.jl")
include("discretemaps/logistic2_unidir.jl")
include("discretemaps/logistic2_bidir.jl")
include("discretemaps/logistic3.jl")
include("discretemaps/logistic4.jl")
include("discretemaps/ulammap.jl")
include("discretemaps/var1.jl")
include("discretemaps/verdes.jl")

include("continuous_systems/chuacircuits_driven.jl")
include("continuous_systems/chuacircuit_nscroll_sine.jl")
include("continuous_systems/hindmarsh_rose.jl")
include("continuous_systems/mediated_link.jl")
include("continuous_systems/lorenz_lorenz_bidir.jl")
include("continuous_systems/lorenz_lorenz_lorenz_bidir_forced.jl")
include("continuous_systems/lorenz_lorenz_lorenz_transitive.jl")
include("continuous_systems/rossler_rossler_bidir.jl")
include("continuous_systems/rossler_rossler_rossler_bidir_forced.jl")
include("continuous_systems/rossler_lorenz.jl")

include("noise.jl")
# We export both the equations of motion and function to generate
# the systems for every example. The equations of motion functions
# are always prepended with `eom_`.
# Some systems (involving multiple lags) do not return instances
# of DiscreteDynamicalSystem, but a Dataset with the finished
# iterated map.

export
###############
# Discrete maps
###############
ar1_unidir,
anishchenko1,
henon2,
henon_triple,
linearmap1,
nonlinear3d,
ulam,
logistic2_unidir,
logistic2_bidir,
logistic3,
logistic4,
nontrivial_pegiun,
var1,
ar1_bidir,
verdes

export
####################
# Continuous systems
###################
chuacircuit_nscroll_sine,
chuacircuits_driven,
hindmarsh_rose,
lorenz_lorenz_lorenz_transitive,
mediated_link,
rossler_rossler_bidir,
rossler_rossler_rossler_bidir_forced,
rossler_lorenz,
lorenz_lorenz_bidir,
lorenz_lorenz_lorenz_bidir_forced

end
