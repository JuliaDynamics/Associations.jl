"""
    AbstractSystemModel{T, N}

An abstract type that represents a set of parameters of type `T` for a dynamical system consisting 
of `N` variables.

Concrete subtypes must implement the following methods:

- **[`get_dt`](@ref)**: Returns the time step.
- **[`get_ui`](@ref)**: Returns the initial condition.
- **[`get_nvars`](@ref)**: Returns the number of variables.

Concrete subtypes *may* implement the following methods:

- **[`rand`](@ref)**. If implemented for a discrete system model of type `SM`, 
    returns an instance of the model with randomised parameters.
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


"""
    rand(::Type{SM}) where {SM <: AbstractSystemModel} -> SM

Generate an instance of a model of type `SM` with randomised parameters.

Parameters are provided as keyword arguments, and must be either as scalars,
distributions or more complicated uncertain values defined as in the 
[UncertainData](https://github.com/kahaaga/UncertainData.jl) package. 

## See also 

Concrete instructions to randomise the different systems are found in the 
following docstrings.

- [`rand(::Type{RosslerLorenzUnidir}`].
"""
function rand(x::AbstractSystemModel) end 

const PT = Union{Number, Distribution, AbstractUncertainValue}