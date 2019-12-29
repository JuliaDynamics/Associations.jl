import StatsBase: mean

"""
    NormalisedPredictiveAsymmetryTest(predictive_test::CausalityTest; f::Number = 1.0)

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
    length(values) == length(ηs) || error("length(values) must equal length(ηs)")
    avg_vals = zeros(Float64, length(ηs[ηs .> 0]))
    
    for (i, η) in enumerate(ηs[ηs .> 0])
        startidx = findfirst(ηs .== -η)
        stopidx  = findlast(ηs .== η)
        avg_vals[i] = mean(values[startidx:stopidx])*f
    end
    
    return avg_vals
end


"""
    normalised_predictive_asymmetry(source, target, test::PredictiveAsymmetryTest{T, N}) where {T <: TransferEntropyCausalityTest, N}

Compute the predictive asymmetry from `source` to `target` using the provided predictive test `p`.
The test can be a [`VisitationFrequencyTest`](@ref), [`TransferOperatorGridTest`](@ref) or 
[`NearestNeighbourMITest`](@ref) and has to be defined for prediction lags symmetrically 
around zero.

## Example 

We'll use the [`logistic4`](@ref) system, which consists of three variables 
and is unidirectionally coupled x -> y -> z.

```julia
# Example data
sys = logistic()
npts = 300
orbit = trajectory(sys, npts, Ttr = 300);
x, y, z = columns(orbit);

# Test setup
η_max = 10
ηs = [-η_max:-1; 1:η_max]
bin = RectangularBinning(floor(Int, npts^(1/4)))
test = VisitationFrequencyTest(ηs = ηs, binning = bin)
pa_test = PredictiveAsymmetryTest(test)

# Analysis
pas_xy = causality(x, y, pa_test)
pas_yx = causality(y, x, pa_test)
pas_yz = causality(y, z, pa_test)
pas_zy = causality(z, y, pa_test)

[pas_xy pas_yx pas_yz pas_zy]
````

If the test correctly picked up the correct directionality, then the values of
of `asymmetries` corresponding to the causal interactions (1st and 3rd column)
should be mostly positive, and those corresponding to non-causal interactions 
should be mostly negative.
"""
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