using Reexport

@reexport module JointDistanceDistribution
    import HypothesisTests: OneSampleTTest, pvalue
    export OneSampleTTest, pvalue 
    include("joint_distance_distribution.jl")
end