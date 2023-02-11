export SystemDefinition, DiscreteDefinition, ContinuousSystem
export system

"""
    SystemDefinition

The abstract type of all system definitions. Abstract subtypes are [`DiscreteDefinition`](@ref)
and [`ContinuousSystem`](@ref). The purpose of concrete implementations is to allow
greater flexibility in the type of systems allowed by `DiscreteDefinition`, and
to streamline randomization of initial conditions and parameters for a particular
set of governing equations.
"""
abstract type SystemDefinition end

"""
    DiscreteDefinition <: SystemDefinition

The supertype of all discrete system definitions.

Why is this type needed? Ideally, an additional `DiscreteDefinition` definition shouldn't
be needed, because we should in principle be able to use `DiscreteDynamicalSystem` directly
for all systems. However, `DiscreteDynamicalSystem` doesn't work
for systems with memory beyond a single time lag. For example, autoregressive systems
of order larger than one are not representable using `DiscreteDynamicalSystem`.

Concrete subtypes of `DiscreteDefinition` are *parameter containers* that are passed
on to [`DiscreteDynamicalSystem`](@ref). They allocate mutable containers that keep
track of past states of state variables that require it. Use [`system`](@ref) to
generate a `DiscreteDynamicalSystem` that can be used with [`trajectory`](@ref).
"""
abstract type DiscreteDefinition <: SystemDefinition end

"""
    ContinuousDefinition <: SystemDefinition

The supertype of all continuous system definitions.
"""
abstract type ContinuousDefinition <: SystemDefinition end

"""
    system(d::DiscreteDefinition)
    system(d::ContinuousDefinition)

Construct a `DiscreteDynamicalSystem` or `ContinuousDynamicalSystem` from `d` that can be
used with trajectory.
"""
function system(d::DiscreteDefinition) end

include("discretemaps/deprecate.jl")
include("discretemaps/AR1Unidir.jl")
include("discretemaps/AR1Bidir.jl")
include("discretemaps/Anishchenko.jl")
include("discretemaps/ChaoticMaps3.jl")
include("discretemaps/henon2.jl")
include("discretemaps/Henon3.jl")
include("discretemaps/Ikeda.jl")
include("discretemaps/LinearMap2.jl")
#include("discretemaps/nonlinear3D_linear_and_nonlinear_coupling.jl")
#include("discretemaps/nontrivial_pegiun.jl")
include("discretemaps/Logistic2Unidir.jl")
include("discretemaps/Logistic2Bidir.jl")
#include("discretemaps/logistic3.jl")
#include("discretemaps/logistic4.jl")
#include("discretemaps/ulammap.jl")
#include("discretemaps/var1.jl")
#include("discretemaps/verdes.jl")

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
