using Reexport

@reexport module SystemModels

"""
    AbstractSystemModel{T, N}

An abstract type that represents a set of parameters of type `T` for a dynamical system consisting 
of `N` variables.

Should implement the following methods:

- **[`get_dt`](@ref)**: Returns the time step.
- **[`get_ui`](@ref)**: Returns the initial condition.
- **[`get_nvars`](@ref)**: Returns the initial condition.
- **[`randomised(::Type{SM}) where SM <: ContiniousSystemModel`](@ref)**. If implemented 
    for a continous system model of type `SM`, returns an instance of the model 
    with randomised parameters.
"""
abstract type AbstractSystemModel{T, N} end 

"""
    DiscreteSystemModel{T, N}

An abstract type that represents a set of model parameters of type T for a discrete dynamical system
consisting of `N` variables.

As a minimum, such types always include a field `dt` (an integer giving the sampling time step, which 
is the `dt` argument to `DynamicalSystems.trajectory`) and a field `ui` field (giving the initial 
condition). Hence, we define a discrete model system by (1) its parameters, (2) its initial condition, 
(3) its sampling rate.

Should implement the following methods:

- **[`get_dt`](@ref)**: Returns the time step.
- **[`get_ui`](@ref)**: Returns the initial condition.
- **[`get_nvars`](@ref)**: Returns the number of variables
- **[`randomised(::Type{SM}) where SM <: DiscreteSystemModel`](@ref)**. If implemented 
    for a discrete system model of type `SM`, returns an instance of the model 
    with randomised parameters.
"""
abstract type DiscreteSystemModel{T, N} <:AbstractSystemModel{T, N} end 

"""
    ContinuousSystemModel{T, N}

An abstract type that represents a set of model parameters of type T for a continous dynamical system,
consisting of `N` variables.

As a minimum, such types always include a field `dt` (giving the sampling time step, which is the 
`dt` argument to `DynamicalSystems.trajectory`) and a field `ui` (giving the initial condition). Hence, 
define a continous model system by (1) its parameters, (2) its initial condition, (3) its sampling rate.

Should implement the following methods:

- **[`get_dt`](@ref)**: Returns the time step.
- **[`get_ui`](@ref)**: Returns the initial condition.
- **[`get_nvars`](@ref)**: Returns the initial condition.
- **[`randomised(::Type{SM}) where SM <: ContiniousSystemModel`](@ref)**. If implemented 
    for a continous system model of type `SM`, returns an instance of the model 
    with randomised parameters.
"""
abstract type ContinuousSystemModel{T, N} <: AbstractSystemModel{T, N} end 

""" Return the integration time step for a system model. """
get_dt(x::AbstractSystemModel) = x.dt

""" Return the initial condition for a system model. """
get_ui(x::AbstractSystemModel) = x.ui

""" Return the number of variables for a system model. """
get_nvars(x::AbstractSystemModel{T, N}) where {T, N} = N

""" 
    randomised(::Type{AbstractSystemModel{T, N}}) where {T, N}

If implemented for the type of system model, return a instance of the model with 
randomised parameters.
"""
function randomised(::Type{SM}) where SM <: ASM where ASM <: AbstractSystemModel
    @error "`randomised` not implemented for system model of type `SM`. Consider defining `randomised(::Type{$(SM)})`"
end

function display_string(x::AbstractSystemModel{T, N}) where {T, N}
    tp = typeof(x)
    fn = fieldnames(typeof(x))
    
    "$tp" * "(" * join(["\n  $(fd) = $(getfield(x, fd))" for fd in fn]) * ")"
end

function Base.show(io::IO, x::AbstractSystemModel{T, N}) where {T, N}
    print(io, display_string(x))
end

export 
AbstractSystemModel, 
ContinuousSystemModel,
DiscreteSystemModel,
randomised,
get_dt,
get_ui,
get_nvars

end