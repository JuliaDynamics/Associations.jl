import StatsBase: mean

"""
    NormalisedPredictiveAsymmetryTest(predictive_test::CausalityTest, f::Number = 1.0)

The parameters for a normalised predictive asymmetry causality test [1]. For the 
non-normalised version, see [`PredictiveAsymmetryTest`](@ref).
    
## Mandatory keywords

- **`predictive_test`**: An instance of a predictive causality test that explicitly 
    uses prediction lags (e.g. [`VisitationFrequencyTest`](@ref) or 
    [`TransferOperatorGridTest`](@ref)). 
- **`f::Number`**: The normalisation factor. 

## About the prediction lags

The prediction lags in the predictive causality test must consist of `n` negative 
integers and `n` positive integers that are symmetric around zero. 

In other words, negative lags  must exactly match the positive lags but with opposite 
sign. The zero lag can be included, but will be ignored, so it is possible to give 
ranges too.

## Examples

Let's define a a transfer entropy test with which to compute transfer entropies
symmetrically around lag zero. We'll use that as input to the predictive asymmetry
test.

```julia 
bin = RectangularBinning(4) # divide each axis into 4 equal-length intervals
ηs = [-5:-1; 1:5] # exclude the zero lag (it is irrelevant for the asymmetry)
test_visitfreq = VisitationFrequencyTest(binning = bin, ηs = ηs)

# Normalising predictive asymmetry values to f = 1.0 times the mean transfer entropy.
NormalisedPredictiveAsymmetryTest(predictive_test = test_visitfreq, f = 1.0)
```

## References 

1. Diego, David, Kristian Agasøster Haaga, Jo Brendryen, and Bjarte Hannisdal. 
    A simple test for causal asymmetry in complex systems. In prep.
"""
Base.@kwdef mutable struct NormalisedPredictiveAsymmetryTest{TEST, N} <: CausalityTest where TEST
    predictive_test::TEST
    f::Number 
    
    # If no threshold is given, use f = 1.0 as default.
    function NormalisedPredictiveAsymmetryTest(test::T) where {T <: TransferEntropyCausalityTest}
        # Check that prediction lags are okay
        verified_prediction_lags(test.ηs)
        N = length(test.ηs[test.ηs .> 0])
        return new{T, N}(test, 1.0)
    end
    
    # If a threshold is given, use it.
    function NormalisedPredictiveAsymmetryTest(test::T; f::Number) where {T <: TransferEntropyCausalityTest}
        # Check that prediction lags are okay
        verified_prediction_lags(test.ηs)
        N = length(test.ηs[test.ηs .> 0])
        return new{T, N}(test, f)
    end
end

""" 
    lagnormalised_statistic(values, ηs, f)

In a cumulative manner, normalise the `values` for the statistic 
to some fraction `f` of their mean over prediction lags `ηs`.
"""
function lagnormalised_statistic(values, ηs, f)
    length(values) == length(ηs) || error("length(vals) must equal length(ηs)")
    avg_vals = zeros(Float64, length(ηs[ηs .> 0]))
    
    for (i, η) in enumerate(ηs[ηs .> 0])
        startidx = findfirst(ηs .== -η)
        stopidx  = findlast(ηs .== η)
        avg_vals[i] = mean(values[startidx:stopidx])*f
    end
    
    return avg_vals
end

function normalised_predictive_asymmetry(source, target, 
        p::NormalisedPredictiveAsymmetryTest{T, N}) where {T <: TransferEntropyCausalityTest, N}

    # Update the test parameters so that we have symmetric prediction lags
    test = update_ηs(p.predictive_test)
    
    # Predictions from source to target 
    preds = causality(source, target, test)
    
    # Normalised predictions
    normalised_preds = lagnormalised_statistic(preds, test.ηs, p.f)
    
    # The number of predictive asymmetries is half the number of prediction lags,
    # which is encoded in the type parameter `N`
    ηs = test.ηs
    As = zeros(Float64, N)

    for (i, η) in enumerate(ηs[ηs .> 0])
        lag_idxs = findfirst(ηs .== -η):findfirst(ηs .== η)
        As[i] = causalbalance(ηs[lag_idxs], preds[lag_idxs]) ./ normalised_preds[i]
    end
    
    return return_predictive_asymmetry(p.predictive_test.ηs, As, N)
end

function normalised_predictive_asymmetry(source, target, cond,
    p::NormalisedPredictiveAsymmetryTest{T, N}) where {T <: TransferEntropyCausalityTest, N}

    # Update the test parameters so that we have symmetric prediction lags
    test = update_ηs(p.predictive_test)

    # Predictions from source to target 
    preds = causality(source, target, cond, test)

    # Normalised predictions
    normalised_preds = lagnormalised_statistic(preds, test.ηs, p.f)

    # The number of predictive asymmetries is half the number of prediction lags,
    # which is encoded in the type parameter `N`
    ηs = test.ηs
    As = zeros(Float64, N)

    for (i, η) in enumerate(ηs[ηs .> 0])
        lag_idxs = findfirst(ηs .== -η):findfirst(ηs .== η)
        As[i] = causalbalance(ηs[lag_idxs], preds[lag_idxs]) ./ normalised_preds[i]
    end

    return return_predictive_asymmetry(p.predictive_test.ηs, As, N)
end


function causality(source::AbstractVector{T}, target::AbstractVector{T}, p::NormalisedPredictiveAsymmetryTest{CT}) where {T<:Real, CT}
    normalised_predictive_asymmetry(source, target, p)
end


function causality(source::AbstractVector{T}, target::AbstractVector{T}, cond::AbstractVector{T}, p::NormalisedPredictiveAsymmetryTest{CT}) where {T<:Real, CT}
    normalised_predictive_asymmetry(source, target, cond, p)
end

export NormalisedPredictiveAsymmetryTest