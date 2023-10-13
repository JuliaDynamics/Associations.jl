import ComplexityMeasures: information 
export information

export MultivariateInformationMeasure
export MultivariateInformationMeasureEstimator
export BivariateInformationMeasure
export BivariateInformationMeasureEstimator

export MutualInformationEstimator
export ConditionalMutualInformationEstimator

# The estimator *always* has the measure definition as the first field with type 
# parameter `M`.
"""
    MultivariateInformationMeasureEstimator

The supertype for all multivariate information measures.

## Concrete subtypes

- [`JointProbababilities`](@ref)
- [`EntropyDecomposition`](@ref)
- [`MIDecomposition`](@ref)

## Abstract subtypes

- [`MutualInformationEstimator`](@ref)
- [`ConditionalMutualInformationEstimator`](@ref)
"""
abstract type MultivariateInformationMeasureEstimator{M} <: InformationMeasureEstimator{M} end
abstract type BivariateInformationMeasureEstimator{M} <: MultivariateInformationMeasureEstimator{M} end

"""
    MutualInformationEstimator

The supertype for dedicated [`MutualInformation`](@ref) estimators.

## Concrete implementations

- [`KSG1`](@ref)
- [`KSG2`](@ref)
- [`GaoOhViswanath`](@ref)
- [`GaoKannanOhViswanath`](@ref)
- [`GaussianMI`](@ref)
"""
abstract type MutualInformationEstimator{M} <: BivariateInformationMeasureEstimator{M} end

"""
    ConditionalMutualInformationEstimator

The supertype for dedicated [`ConditionalMutualInformation`](@ref) estimators.

## Concrete implementations

- [`FPVP`](@ref)
- [`GaussianCMI`](@ref)
- [`MesnerShalizi`](@ref)
- [`Rahimzamani`](@ref)
- [`PoczosSchneiderCMI`](@ref)
"""
abstract type ConditionalMutualInformationEstimator{M} <: MultivariateInformationMeasureEstimator{M} end

"""
    MultivariateInformationMeasure <: AssociationMeasure

The supertype for all multivariate information-based measure definitions.

## Implementations

- [`JointEntropy`](@ref) (its concrete subtypes)
- [`ConditionalEntropy`](@ref) (its concrete subtypes)
- [`MutualInformation`](@ref) (its concrete subtypes)
- [`ConditionalMutualInformation`](@ref) (its concrete subtypes)
- [`TransferEntropy`](@ref) (its concrete subtypes)
- [`PartialMutualInformation`](@ref)
"""
abstract type MultivariateInformationMeasure <: AssociationMeasure end

"""
    BivariateInformationMeasure <: MultivariateInformationMeasure

The supertype of all bivariate information measure definitions.
"""
abstract type BivariateInformationMeasure <: MultivariateInformationMeasure end

min_inputs_vars(::BivariateInformationMeasure) = 2
max_inputs_vars(::BivariateInformationMeasure) = 2

function verify_number_of_inputs_vars(est::MultivariateInformationMeasureEstimator, n)
    return verify_number_of_inputs_vars(est.definition, n)
end

# For internal use only
function estimate(est::MultivariateInformationMeasureEstimator, x...)
    return information(est, x...)
end