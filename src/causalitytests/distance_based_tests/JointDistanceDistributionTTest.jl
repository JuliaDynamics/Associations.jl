



"""
JointDistanceDistributionTTest(; distance_metric = SqEuclidean(), B::Int = 10,
    D::Int = 2, τ::Int = 1, 
    hypothesis_test::OneSampleTTest = OneSampleTTest,
    μ0 = 0.0)

The parameters for a joint distance distribution [1] analysis.

## Optional keyword arguments 

- **`distance_metric`**: The distance metric used to compute distances. Has to be a instance of a 
valid distance metric from `Distances.jl`. Defaults to `SqEuclidean()`.

- **`B::Int`**: The number of equidistant subintervals to divide the interval `[0, 1]` into
when comparing the normalised distances. 

- **`D::Int`**: The dimension of the delay reconstructions.

- **`τ::Int`**: The delay of the delay reconstructions.

- **`μ0`**: The hypothetical mean value of the joint distance distribution if there 
is no coupling between `x` and `y` (default is `μ0 = 0.0`).

- **`hypothesis_test`**: A `OneSampleTTest` to test whether the joint distance distribution 
is skewed towards positive values.

## References 

[1] Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from multivariate flows by 
the joint distance distribution." Chaos: An Interdisciplinary Journal of Nonlinear 
Science 28.7 (2018): 075302.
"""
Base.@kwdef struct JointDistanceDistributionTTest <: JointDistancesCausalityTest
    """ The distance metric. """
    distance_metric = SqEuclidean() 

    """ The number of subintervals. """
    B::Int = 10

    """ The dimension of the delay reconstructions. """
    D::Int = 2

    """ The delay reconstruction lag. """
    τ::Int = 1

    """ The hypothesis test. Initialize with `nothing` if no test is to be applied."""
    hypothesis_test::Type{OneSampleTTest} = OneSampleTTest

    """ The default `μ0` value for the default `OneSampleTTest` hypothesis test. """
    μ0 = 0.0
end

function causality(source, target, p::JointDistanceDistributionTTest)
    joint_distance_distribution(p.hypothesis_test, source, target,
        p.distance_metric, 
        p.B,
        p.D,
        p.τ,
        p.μ0)
end

################################################################
# Integration with UncertainData.jl
################################################################
function causality(source::Vector{<:UVT1}, target::Vector{<:UVT2}, p::JointDistanceDistributionTTest) where {
        UVT1<:AbstractUncertainValue, 
        UVT2<:AbstractUncertainValue}
    causality(resample.(source), resample.(target), p)
end

function causality(source, target::Vector{<:UVT}, p::JointDistanceDistributionTTest) where { 
        UVT<:AbstractUncertainValue}
    causality(source, resample.(target), p)
end

function causality(source::Vector{<:UVT}, target, p::JointDistanceDistributionTTest) where {
        UVT<:AbstractUncertainValue}
    causality(resample.(source), target, p)
end

function causality(source::UVT1, target::UVT2, p::JointDistanceDistributionTTest) where {
        UVT1<:AbstractUncertainValueDataset, 
        UVT2<:AbstractUncertainValueDataset}
    causality(resample.(source), resample.(target), p)
end

function causality(source, target::UVT, p::JointDistanceDistributionTTest) where { 
        UVT<:AbstractUncertainValueDataset}
    causality(source, resample.(target), p)
end

function causality(source::UVT, target, p::JointDistanceDistributionTTest) where {
        UVT<:AbstractUncertainValueDataset}
    causality(resample.(source), target, p)
end

export JointDistanceDistributionTTest