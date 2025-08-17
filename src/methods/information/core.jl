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

- [`JointProbabilities`](@ref)
- [`EntropyDecomposition`](@ref)
- [`MIDecomposition`](@ref)
- [`CMIDecomposition`](@ref)

## Dedicated implementations

[`MutualInformationEstimator`](@ref)s:

- [`KraskovStögbauerGrassberger2`](@ref)
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

- [`KraskovStögbauerGrassberger2`](@ref)
- [`KraskovStögbauerGrassberger2`](@ref)
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

## Definition

Following [Datseris2024](@citet), we define a multivariate information measure as *any functional 
of a multidimensional probability mass functions (PMFs) or multidimensional probability density*.

## Implementations

[`JointEntropy`](@ref) definitions:

- [`JointEntropyShannon`](@ref)
- [`JointEntropyRenyi`](@ref)
- [`JointEntropyTsallis`](@ref)

[`ConditionalEntropy`](@ref) definitions:

- [`ConditionalEntropyShannon`](@ref)
- [`ConditionalEntropyTsallisAbe`](@ref)
- [`ConditionalEntropyTsallisFuruichi`](@ref)

[`DivergenceOrDistance`](@ref) definitions:

- [`HellingerDistance`](@ref)
- [`KLDivergence`](@ref)
- [`RenyiDivergence`](@ref)
- [`VariationDistance`](@ref)

[`MutualInformation`](@ref) definitions:

- [`MIShannon`](@ref)
- [`MIRenyiJizba`](@ref)
- [`MIRenyiSarbu`](@ref)
- [`MITsallisMartin`](@ref)
- [`MITsallisFuruichi`](@ref)

[`ConditionalMutualInformation`](@ref) definitions:

- [`CMIShannon`](@ref)
- [`CMITsallisPapapetrou`](@ref)
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

function size_match(measure::BivariateInformationMeasure, px::Probabilities, py::Probabilities)
    size(px) == size(py) || throw(DimensionMismatch("px and py must have the same size"))
end

"""
    information(est::MultivariateInformationMeasureEstimator, x...)

Estimate some [`MultivariateInformationMeasure`](@ref) on input data `x...`,
using the given [`MultivariateInformationMeasureEstimator`](@ref).

This is just a convenience wrapper around [`association`](@ref)`(est, x...)`.
"""
function information(est::MultivariateInformationMeasureEstimator, x...)
    return association(est, x...)
end