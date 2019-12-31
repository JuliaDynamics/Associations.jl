function get_return_type(test::TransferEntropyCausalityTest)
    T = typeof(1.0)

    if test.ηs isa Number 
        return T
    else test.ηs isa Vector
        return Vector{T}
    end
end

function get_return_type(test::CrossMappingTest)
    T = typeof(1.0)

    return Vector{T}
end

function get_return_type(test::ConvergentCrossMappingTest)
    T = typeof(1.0)

    return Vector{Vector{T}}
end

function get_return_type(test::JointDistanceDistributionTest)
    T = typeof(1.0)

    return Vector{T}
end

function get_return_type(test::JointDistanceDistributionTTest)
    return OneSampleTTest
end

function get_return_type(test::SMeasureTest)
    T = typeof(1.0)
    return Vector{T}
end


function get_return_type(test::AbstractPredictiveAsymmetryTest)
    T = typeof(1.0)
    return Vector{T}
end

export get_return_type