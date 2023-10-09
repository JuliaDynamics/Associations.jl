using Accessors: @set

export estimator_with_overridden_parameters

# For internal use only.
"""
    estimator_with_overridden_parameters(definition, lower_level_estimator) → e::typeof(lower_level_estimator)

Given some higher-level `definition` of an information measure, which is to be 
estimated using some `lower_level_estimator`, return a modified version of 
the estimator in which its parameter have been overriden by any overlapping
parameters from the `defintiion`.

This method is explicitly extended for each possible decomposition.
"""
function estimator_with_overridden_parameters(definition, lower_level_estimator) end

const TSALLIS_MULTIVARIATE_MEASURES = Union{
    CMITsallis, 
    MITsallisFuruichi, MITsallisMartin,
    CETsallisAbe, CETsallisFuruichi,
    JointEntropyTsallis,
}

function estimator_with_overridden_parameters(defintion::CMITsallis, 
        est::InformationMeasureEstimator{<:Tsallis})
    modified_est = @set est.base = definition.base
    modified_est = @set est.q = definition.q

    return modified_est
end