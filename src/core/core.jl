using DelayEmbeddings: AbstractStateSpaceSet
using ComplexityMeasures: ProbabilitiesEstimator
const VectorOrStateSpaceSet{D, T} = Union{AbstractVector{T}, AbstractStateSpaceSet{D, T}} where {D, T}
const ArrayOrStateSpaceSet{D, T, N} = Union{AbstractArray{T, N}, AbstractStateSpaceSet{D, T}} where {D, T, N}
const ENTROPY_ESTS = Union{DifferentialInfoEstimator, DiscreteInfoEstimator}

export AssociationMeasure
export DirectedAssociationMeasure
export Discretization
# The type parameter `N` indicates the number of input datasets to be discretized.
"""
    Discretization

The supertype of all discretization schemes.

## Description 

A fundamental operation when computing multivariate information measures from data 
is *discretization*. There are many ways of discretizing multiple input datasets. We 
offer two main ways of doing so.

## Concrete implementations

- [`CodifyVariables`](@ref)
- [`CodifyPoints`](@ref)
"""
abstract type Discretization{N} end

# Any non-bivariate association measures must implement:
# - [`min_inputs_vars`](@ref).
# - [`max_inputs_vars`](@ref).
"""
    AssociationMeasure

The supertype of all association measures.
"""
abstract type AssociationMeasure end

abstract type DirectedAssociationMeasure <: AssociationMeasure end

# Just use ComplexityMeasures.convert_logunit when it is released.
"""
    _convert_logunit(h_a::Real, , to) → h_b

Convert a number `h_a` computed with logarithms to base `a` to an entropy `h_b` computed
with logarithms to base `b`. This can be used to convert the "unit" of an entropy.
"""
function _convert_logunit(h::Real, base_from, base_to)
    h / log(base_from, base_to)
end

# Default to bivariate measures. Other measures override it.
"""
    min_inputs_vars(m::AssociationMeasure) → nmin::Int

Return the minimum number of variables is that the measure can be computed for.

For example, [`CMIShannon`](@ref) requires 3 input variables.
"""
min_inputs_vars(m::AssociationMeasure) = 2

# Default to bivariate measures. Other measures override it.

"""
    max_inputs_vars(m::AssociationMeasure) → nmax::Int

Return the maximum number of variables is that the measure can be computed for.

For example, [`MIShannon`](@ref) cannot be computed for more than 2 variables.
"""
max_inputs_vars(m::AssociationMeasure) = 2

function verify_number_of_inputs_vars(measure::AssociationMeasure, n::Int)
    T = typeof(measure)
    nmin = min_inputs_vars(measure)
    if n < nmin
        throw(ArgumentError("$T requires at least $nmin inputs. Got $n inputs."))
    end

    nmax = max_inputs_vars(measure)
    if n > nmax
        throw(ArgumentError("$T accepts a maximum of $nmax inputs. Got $n inputs."))
    end
end

# Concrete ways of encoding multivariate data, each defined as a type.
include("encoding/codify_points.jl")
include("encoding/codify_variables.jl")

# Counting and probabilities (contingency tables and probabilities for multivariate data)
include("counts/counts.jl")
include("probabilities/probabilities.jl")

# Estimating counts/probabilities from multiple input datasets.
include("estimation_per_point.jl")
include("estimation_per_variable.jl")
