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

const RENYI_MULTIVARIATE_MEASURES = Union{
    CMIRenyiPoczos, CMIRenyiSarbu, CMIRenyiJizba,
    MIRenyiJizba, MIRenyiSarbu,
    JointEntropyRenyi,
}

const SHANNON_MULTIVARIATE_MEASURES = Union{
    CMIShannon,
    MIShannon,
    CEShannon,
    JointEntropyShannon,
}


function estimator_with_overridden_parameters(
        definition::TSALLIS_MULTIVARIATE_MEASURES, 
        est::InformationMeasureEstimator{<:Tsallis}
    )
    # The low-level definition
    lowdef = est.definition
   
    # Update the low-level definition. Have to do this step-wise. Ugly, but works.
    modified_lowdef = @set lowdef.base = definition.base # update `base` field
    modified_lowdef = @set modified_lowdef.q = definition.q # update `q` field

    # Set the definition for the low-level estimator to the updated definition.
    modified_est = @set est.definition = modified_lowdef
    
    return modified_est
end

function estimator_with_overridden_parameters(
        definition::RENYI_MULTIVARIATE_MEASURES, 
        est::InformationMeasureEstimator{<:Renyi}
    )
    lowdef = est.definition
    modified_lowdef = @set lowdef.base = definition.base # update `base` field
    modified_lowdef = @set modified_lowdef.q = definition.q # update `q` field
    modified_est = @set est.definition = modified_lowdef
    return modified_est
end

function estimator_with_overridden_parameters(
        definition::SHANNON_MULTIVARIATE_MEASURES, 
        est::InformationMeasureEstimator{<:Shannon}
    )
    lowdef = est.definition
    modified_lowdef = @set lowdef.base = definition.base # update `base` field
    modified_est = @set est.definition = modified_lowdef
    return modified_est
end