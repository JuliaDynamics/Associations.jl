"""
    AbstractSystemModel{T, N}

An abstract type that represents a set of parameters of type `T` for a dynamical system consisting 
of `N` variables.

Concrete subtypes should implement the following methods:

- **[`get_dt`](@ref)**: Returns the time step.
- **[`get_ui`](@ref)**: Returns the initial condition.
- **[`get_nvars`](@ref)**: Returns the initial condition.
- **[`randomised(::Type{SM}) where SM <: ContiniousSystemModel`](@ref)**. If implemented 
    for a continous system model of type `SM`, returns an instance of the model 
    with randomised parameters.
"""
abstract type AbstractSystemModel{T, N} end 


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


export get_dt, get_ui, get_nvars, randomised