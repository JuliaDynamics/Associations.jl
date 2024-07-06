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

Concrete subtypes are given as input to [`association`](@ref).

## Concrete implementations 

- [`PearsonCorrelation`](@ref)
- [`PartialCorrelation`](@ref)
- [`DistanceCorrelation`](@ref)
- [`SMeasure`](@ref)
- [`LMeasure`](@ref)
- [`HMeasure`](@ref)
- [`MMeasure`](@ref)
- [`JointDistanceDistribution`](@ref).

## Abstract implementations

The docstrings for the abstract types below list concrete implementations.

- [`MultivariateInformationMeasure`](@ref)
- [`CrossmapMeasure`](@ref)
"""
abstract type AssociationMeasure end

"""
    AssociationMeasureEstimator

The supertype of all association measure estimators.

Concrete subtypes are given as input to [`association`](@ref).

## Abstract subtypes

See the documentation of the abstract types below for concrete implementations.

- [`MultivariateInformationMeasureEstimator`](@ref)
- [`CrossmapEstimator`](@ref)
"""
abstract type AssociationMeasureEstimator end

"""
    association(estimator::AssociationMeasureEstimator, x, y, [z, ...]) → r
    association(definition::AssociationMeasure, x, y, [z, ...]) → r

Estimate the (conditional) association between input variables `x, y, z, …` using 
the given `estimator` or `definition`.  

The type of the return value `r` depends on the `measure`/`estimator`.

## Compatible definitions and estimators

Some [`AssociationMeasure`](@ref)s have no estimators and are given directly.
For other association measures, you need to provide an [`AssociationMeasureEstimator`](@ref)
 of some kind, which takes the definition as its first argument.

## Examples

```julia
using CausalityTools
x, y, z = rand(1000), rand(1000), rand(1000)

# Pairwise association using different measures
association(DistanceCorrelation(), x, y)
association(PartialCorrelation(), x, y)
association(JointProbabilities(ConditionalEntropyTsallisAbe(), ValueBinning(3)), x, y)
association(JointProbabilities(JointEntropyShannon(), Dispersion(c = 3, m = 2)), x, y)
association(EntropyDecomposition(MIShannon(), PlugIn(Shannon()), OrdinalPatterns(m=3)), x, y)
association(KSG2(MIShannon(base = 2)), x, y)

# Conditional association using different measures
association(JointProbabilities(PartialMutualInformation(), OrdinalPatterns(m=3)), x, y, z)
association(FPVP(CMIShannon(base = 2)), x, y, z)
```
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