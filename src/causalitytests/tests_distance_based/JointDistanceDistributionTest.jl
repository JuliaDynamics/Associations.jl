
"""
    JointDistanceDistributionTest(; distance_metric = SqEuclidean(), B::Int = 10,
        D::Int = 2, τ::Int = 1)

The parameters for a joint distance distribution [1] analysis.

## Optional keyword arguments 

- **`distance_metric`**: The distance metric used to compute distances. Has to be a instance of a 
    valid distance metric from `Distances.jl`. Defaults to `SqEuclidean()`.

- **`B::Int`**: The number of equidistant subintervals to divide the interval `[0, 1]` into
    when comparing the normalised distances. 

- **`D::Int`**: The dimension of the delay reconstructions.

- **`τ::Int`**: The delay of the delay reconstructions.

## References 

[1] Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from multivariate flows by 
the joint distance distribution." Chaos: An Interdisciplinary Journal of Nonlinear 
Science 28.7 (2018): 075302.

"""
Base.@kwdef struct JointDistanceDistributionTest{N, M} <: JointDistancesCausalityTest{N} where M
    """ The distance metric. """
    distance_metric::M = SqEuclidean() 
    
    """ The number of subintervals. """
    B::Int = 10
    
    """ The dimension of the delay reconstructions. """
    D::Int = 2
    
    """ The delay reconstruction lag. """
    τ::Int = 1

    function JointDistanceDistributionTest(distance_metric::M, B::Int, D::Int, τ::Int) where M
        # The number of returned elements. The test itself returns a distribution of 
        # `B` different numbers, but we're applying the hypothesis test to that, so 
        # only a single element is returned.
        N = B
        new{N, M}(distance_metric, B, D, τ)
    end
end


function causality(source::AbstractVector{T}, target::AbstractVector{T}, p::JointDistanceDistributionTest) where {T <: Real}
    joint_distance_distribution(source, target;
        distance_metric = p.distance_metric, 
        B = p.B, 
        D = p.D, 
        τ = p.τ)
end

export JointDistanceDistributionTest