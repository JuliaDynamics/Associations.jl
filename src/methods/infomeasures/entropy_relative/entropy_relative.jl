export RelativeEntropyEstimator
export entropy_relative

abstract type RelativeEntropyEstimator end

include("tsallis/entropy_relative_tsallis.jl")


function entropy_relative(e::Entropy, est::RelativeEntropyEstimator,
        x::Vector_or_Dataset, y::Vector_or_Dataset)
    throw(ArgumentError(
        "Relative $typeof(e) entropy not defined for $typeof(est) estimator"
    ))
end


include("estimators/Wang.jl")
