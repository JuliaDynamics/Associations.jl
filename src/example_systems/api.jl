export SystemDefinition, DiscreteDefinition, ContinuousDefinition, LaggedDiscreteDefinition
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
"""
abstract type DiscreteDefinition <: SystemDefinition end

"""
    ContinuousDefinition <: SystemDefinition

The supertype of all continuous system definitions.
"""
abstract type ContinuousDefinition <: SystemDefinition end

"""
    system(definition::DiscreteDefinition) → s::DiscreteDynamicalSystem
    system(definition::ContinuousDefinition) → s::ContinuousDynamicalSystem

Initialize a dynamical system from `definition`.
"""
function system(d::SystemDefinition) end

################################################################
# Internal type for lagged-system-specific dispatch.
################################################################
"""
    LaggedDiscreteDefinition <: SystemDefinition

The supertype of definitions for discrete systems with lag larger than 1.

Why is this type needed? Ideally, an additional definition shouldn't
be needed, because we should in principle be able to use `DiscreteDynamicalSystem` directly
for all systems. However, `DiscreteDynamicalSystem` doesn't work
for systems with memory beyond a single time lag. For example, autoregressive systems
of order larger than one are not representable using `DiscreteDynamicalSystem`.

Concrete subtypes of `DiscreteDefinition` are *parameter containers* that are passed
on to [`DiscreteDynamicalSystem`](@ref). They allocate mutable containers that keep
track of past states of state variables that require it. Use [`system`](@ref) to
generate a `DiscreteDynamicalSystem` that can be used with [`trajectory`](@ref).

## Implementation details

Concrete implementations must fulfill the below criteria.

- Subtypes must implement a `past_states` field, which is
    a `SVector{N, MVector{L, Int}}`, where `N` is the number of variables. For type stability,
    `L` states are tracked for *all* variables, even though the maximum lag may only occur
    for one of the variables.
- The first type parameter of subtypes must be `P`, which keeps track of the type of
    `past_states`.

For an example, see the source code for [`Peguin2`](@ref).
"""
abstract type LaggedDiscreteDefinition{P} <: SystemDefinition end

"""
    update_states!(s::LaggedDiscreteDefinition, xnew::SVector{D})

Given `xnew` (the new current state of a system), update the past states of `s`.
"""
function update_states!(def::LaggedDiscreteDefinition{SVector{N, MVector{D, T}}},
        xnew::SVector{N}) where {N, D, T}
    D >= 2 || error("Memory vector for LaggedDiscreteDefinition must have length at least 2.")
    for var in 1:N
        # A LaggedDiscreteDefinition always has a memory vector of length at least 2.
        # Otherwise, it should be a regular DiscreteDefinition.
        for k in D:-1:2
            def.past_states[var][k] = def.past_states[var][k - 1]
        end
        def.past_states[var][1] = xnew[var]
    end
end
