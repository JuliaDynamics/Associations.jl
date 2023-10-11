export MultivariateInformationMeasure
export MultivariateInformationMeasureEstimator
export BivariateInformationMeasure
export BivariateInformationMeasureEstimator

"""
    MultivariateInformationMeasure <: AssociationMeasure

The supertype for all multivariate information-based measure definitions.
"""
abstract type MultivariateInformationMeasure <: AssociationMeasure end

"""
    BivariateInformationMeasure <: AssociationMeasure

The supertype of all bivariate information measures (defined here as functionals
of probability mass functions or probability densities).

## Implements

- [`information`](@ref). Used to compute the measure given given relevant input
    probabilities.

## Concrete implementations

### Discrepancy/divergence/closeness measures

A multitude of measures exist to quantify the discrepancy/divergence/closeness between two
probability distributions. Some of these are proper metrics, others are not.
They have in common that they aim to quantify how far from each other two
probabilities distributions are.

- [`KLDivergence`](@ref)
- [`RenyiDivergence`](@ref)
- [`HellingerDistance`](@ref)

The supertype for all bivariate information-based measure definitions.
"""
abstract type BivariateInformationMeasure <: MultivariateInformationMeasure end

min_inputs_vars(::BivariateInformationMeasure) = 2
max_inputs_vars(::BivariateInformationMeasure) = 2

abstract type MultivariateInformationMeasureEstimator end
abstract type BivariateInformationMeasureEstimator end
