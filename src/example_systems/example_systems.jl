export SystemDefinition, DiscreteSystem, ContinuousSystem
export system

"""
    SystemDefinition

The abstract type of all system definitions. Abstract subtypes are [`DiscreteSystem`](@ref)
and [`ContinuousSystem`](@ref). The purpose of concrete implementations is to allow
greater flexibility in the type of systems allowed by `DiscreteSystem`, and
to streamline randomization of initial conditions and parameters for a particular
set of governing equations.
"""
abstract type SystemDefinition end

"""
    DiscreteSystem <: SystemDefinition

The supertype of all discrete system definitions.

Why is this type needed? Ideally, an additional `DiscreteSystem` definition shouldn't
be needed, because we should in principle be able to use `DiscreteDynamicalSystem` directly
for all systems. However, `DiscreteDynamicalSystem` doesn't work
for systems with memory beyond a single time lag. For example, autoregressive systems
of order larger than one are not representable using `DiscreteDynamicalSystem`.

Concrete subtypes of `DiscreteSystem` are *parameter containers* that are passed
on to [`DiscreteDynamicalSystem`](@ref). They allocate mutable containers that keep
track of past states of state variables that require it. Use [`system`](@ref) to
generate a `DiscreteDynamicalSystem` that can be used with [`trajectory`](@ref).
"""
abstract type DiscreteSystem end

"""
    ContinuousSystem <: SystemDefinition

The supertype of all continuous system definitions.
"""
abstract type ContinuousSystem end

"""
    system(d::DiscreteSystem)
    system(d::ContinuousSystem)

Construct a `DiscreteDynamicalSystem` or `ContinuousDynamicalSystem` from `d` that can be
used with trajectory.
"""
function system(d::DiscreteSystem) end

include("discretemaps/deprecate.jl")
include("discretemaps/ar1.jl")
include("discretemaps/ar1_bidir.jl")
include("discretemaps/anishchenko1.jl")
include("discretemaps/henon2.jl")
include("discretemaps/Henon3.jl")
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
include("continuous_systems/repressilator.jl")
include("continuous_systems/rossler_rossler_bidir.jl")
include("continuous_systems/rossler_rossler_rossler_bidir_forced.jl")
include("continuous_systems/rossler_lorenz.jl")

include("noise.jl")
