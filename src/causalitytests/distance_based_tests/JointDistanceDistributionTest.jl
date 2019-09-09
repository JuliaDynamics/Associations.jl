import Distances: SqEuclidean
import HypothesisTests: OneSampleTTest

import Distances: SqEuclidean
import HypothesisTests: OneSampleTTest

"""
    JointDistancesTest

The supertype of all joint distance distribution tests.
"""
abstract type JointDistancesTest <: DistanceBasedCausalityTest end

"""
    JointDistanceDistributionTest(; distance_metric = SqEuclidean(), B::Int = 10,
        D::Int = 2, τ::Int = 1)

The parameters for a joint distance distribution analysis.
"""
Base.@kwdef struct JointDistanceDistributionTest <: JointDistancesTest
    """ The distance metric. """
    distance_metric = SqEuclidean() 
    
    """ The number of subintervals. """
    B::Int = 10
    
    """ The dimension of the delay reconstructions. """
    D::Int = 2
    
    """ The delay reconstruction lag. """
    τ::Int = 1
end

"""
    JointDistanceDistributionTTest(; distance_metric = SqEuclidean(), B::Int = 10,
        D::Int = 2, τ::Int = 1, 
        hypothesis_test::OneSampleTTest = OneSampleTTest,
        μ0 = 0.0)

The parameters for a joint distance distribution analysis.
"""
Base.@kwdef struct JointDistanceDistributionTTest <: JointDistancesTest
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

function causality(source, target, p::JointDistanceDistributionTest)
    joint_distance_distribution(source, target;
        distance_metric = p.distance_metric, 
        B = p.B, 
        D = p.D, 
        τ = p.τ)
end


function causality(source, target, p::JointDistanceDistributionTTest)
    joint_distance_distribution(p.hypothesis_test, source, target,
        p.distance_metric, p.B, p.D, p.τ, p.μ0)
end

export 
    JointDistanceDistributionTest, 
    JointDistanceDistributionTTest