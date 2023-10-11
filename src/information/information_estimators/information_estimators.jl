include("Contingency.jl")
include("DiscreteDecomposition.jl")
include("DifferentialDecomposition.jl")
include("EntropyDecomposition.jl")
abstract type MutualInformationEstimator{M} <: BivariateInformationMeasureEstimator end