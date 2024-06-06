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

The supertype for all estimators of multivariate information measures.

## Generic implementations

- [`JointProbababilities`](@ref)
- [`EntropyDecomposition`](@ref)
- [`MIDecomposition`](@ref)
- [`CMIDecomposition`](@ref)

## Dedicated implementations

[`MutualInformationEstimator`](@ref)s:

- [`KraskovStögbauerGrassberger1`](@ref)
- [`KraskovStögbauerGrassberger2`](@ref)
- [`GaoOhViswanath`](@ref)
- [`GaoKannanOhViswanath`](@ref)
- [`GaussianMI`](@ref)

[`ConditionalMutualInformationEstimator`](@ref)s:

- [`FPVP`](@ref)
- [`MesnerShalizi`](@ref)
- [`Rahimzamani`](@ref)
- [`PoczosSchneiderCMI`](@ref)
- [`GaussianCMI`](@ref)

[`TransferEntropyEstimator`](@ref)s:

- [`Zhu1`](@ref)
- [`Lindner`](@ref)
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

[`JointEntropy`](@ref) definitions:

- [`JointEntropyShannon`](@ref)
- [`JointEntropyRenyi`](@ref)
- [`JointEntropyTsallis`](@ref)

[`ConditionalEntropy`](@ref) definitions:

- [`ConditionalEntropyShannon`](@ref)
- [`ConditionalEntropyShannon`](@ref)

[`DivergenceOrDistance`](@ref) definitions:

- [`HellingerDistance`](@ref)
- [`KLDivergence`](@ref)
- [`RenyiDivergence`](@ref)
- [`VariationDistance`](@ref)

[`MutualInformation`](@ref) definitions:

- [`MIShannon`](@ref)
- [`MIRenyiJizba`](@ref)
- [`MIRenyiMartin`](@ref)
- [`MITsallisAbe`](@ref)
- [`MITsallisFuruchi`](@ref)

[`ConditionalMutualInformation`](@ref) definitions:

- [`CMIShannon`](@ref)
- [`CMITsallis`](@ref)
- [`CMIRenyiJizba`](@ref)
- [`CMIRenyiPoczos`](@ref)
- [`CMIRenyiSarbu`](@ref)

[`TransferEntropy`](@ref) definitions:

- [`TEShannon`](@ref)
- [`TERenyiJizba`](@ref)

Other definitions:

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