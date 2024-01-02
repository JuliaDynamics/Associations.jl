using ComplexityMeasures: DifferentialInfoEstimator, DiscreteInfoEstimator
using StateSpaceSets: AbstractStateSpaceSet

export AssociationMeasure
export AssociationMeasureEstimator
export association

const VectorOrStateSpaceSet{D, T} = Union{AbstractVector{T}, AbstractStateSpaceSet{D, T}} where {D, T}
const ArrayOrStateSpaceSet{D, T, N} = Union{AbstractArray{T, N}, AbstractStateSpaceSet{D, T}} where {D, T, N}
const ENTROPY_ESTS = Union{DifferentialInfoEstimator, DiscreteInfoEstimator}

"""
    AssociationMeasure

The supertype of all association measures.
"""
abstract type AssociationMeasure end

"""
    AssociationMeasureEstimator

The supertype of all association measures.
"""
abstract type AssociationMeasureEstimator end

"""
    association(definition::AssociationMeasure, x, y, [z, ...])
    association(est::AssociationMeasureEstimator, x, y, [z, ...])

Compute the association between input variables `x, y, z, …` according to the 
estimator `est` (or directly from the `definition` if possible).

The return type depends on the definition/estimator.

## Compatible definitions

Some measures are estimated directly (they have no estimator types). These can be 
given to `association` directly. They are:

- [`PearsonCorrelation`](@ref)
- [`PartialCorrelation`](@ref)
- [`DistanceCorrelation`](@ref)
- [`SMeasure`](@ref)
- [`LMeasure`](@ref)
- [`HMeasure`](@ref)
- [`MMeasure`](@ref)
- [`JointDistanceDistribution`](@ref).

## Compatible estimators

For all other association measures, you need to provide an
[`AssociationMeasureEstimator`](@ref) of some kind. It takes the definition as its 
first argument.

### Cross mappings

- [`RandomVectors`](@ref), [`RandomSegments`](@ref), [`ExpandingSegments`](@ref), which 
    all take either a [`ConvergentCrossMapping`](@ref) or
    [`PairwiseAsymmetricEmbedding`](@ref) definition as their first arg


"""
function association(est, x...) end

# Just use ComplexityMeasures.convert_logunit when it is released.
"""
    _convert_logunit(h_a::Real, , to) → h_b

Convert a number `h_a` computed with logarithms to base `a` to an entropy `h_b` computed
with logarithms to base `b`. This can be used to convert the "unit" of e.g. an entropy.
"""
function _convert_logunit(h::Real, base_from, base_to)
    h / log(base_from, base_to)
end

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