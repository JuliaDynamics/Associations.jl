include("Contingency.jl")
include("DiscreteDecomposition.jl")
include("DifferentialDecomposition.jl")

include("EntropyDecomposition.jl")
include("MIDecomposition.jl")

abstract type MutualInformationEstimator{M} <: BivariateInformationMeasureEstimator end
abstract type ConditionalMutualInformationEstimator{M} <: MultivariateInformationMeasureEstimator end
