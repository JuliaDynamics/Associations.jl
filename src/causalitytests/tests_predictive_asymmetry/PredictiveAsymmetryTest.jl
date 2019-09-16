

function verified_prediction_lags(lags)
    ηs = sort(lags[lags .!= 0])

    if length(ηs) % 2 != 0
        throw(ArgumentError("Need an even number of lags (as many negative as positive prediction lags)"))
    end

    ηs_neg = ηs[ηs .< 0]
    ηs_pos = ηs[ηs .> 0]

    if length(ηs_neg) != length(ηs_pos) 
        throw(ArgumentError("There must be as many negative as positive prediction lags. Got $ηs_neg and $ηs_pos"))
    end

    ηs_neg_sort = sort(abs.(ηs_neg))
    ηs_pos_sort = sort(ηs_pos)
    if !(all(ηs_neg_sort .== ηs_pos_sort))
        throw(ArgumentError("Negative amd positive prediction lags must be symmetric. Got $ηs_neg and $ηs_pos"))
    end 

    return ηs
end

function verified_prediction_lags(lag::Int)
    # If only one η is provided (either positive or negative), just make a prediction 
    # for that lag.
    ηs = [-abs(lag), abs(lag)]
end

"""
    PredictiveAsymmetryTest(predictive_test::CausalityTest)

The parameters for a predictive asymmetry causality test [1]. 
    
## Mandatory keywords

- **`predictive_test`**: An instance of a predictive causality test that explicitly 
    uses prediction lags (e.g. [`VisitationFrequencyTest`](@ref) or 
    [`TransferOperatorGridTest`](@ref)). 

## About the prediction lags

The prediction lags in the predictive causality test must consist of `n` negative 
integers and `n` positive integers that are symmetric around zero. 

In other words, negative lags  must exactly match the positive lags but with opposite 
sign. The zero lag can be included, but will be ignored, so it is possible to give 
ranges too.

## Examples

```julia 
test_visitfreq = VisitationFrequencyTest(ηs = [-5, -4, -2, -1, 0, 1, 2, 4, 5])
test_transferoperator = TransferOperatorGridTest(ηs = -3:3)

# Note that `predictive_test` is a *mandatory* keyword.
PredictiveAsymmetryTest(predictive_test = test_visitfreq)
PredictiveAsymmetryTest(predictive_test = test_transferoperator)
```

## References 

1. Diego, David, Kristian Agasøster Haaga, Jo Brendryen, and Bjarte Hannisdal. 
    A simple test for causal asymmetry in complex systems. In prep.
"""
Base.@kwdef struct PredictiveAsymmetryTest{TEST} <: CausalityTest where TEST
    predictive_test::TEST

    function PredictiveAsymmetryTest(test::T) where {T <: TransferEntropyCausalityTest}
        # Check that prediction lags are okay
        verified_prediction_lags(test.ηs)
        return new{T}(test)
    end
end

function causalbalance(lags, tes)
    posinds = findall(lags .> 0)
    neginds = findall(lags .< 0)
    sum(tes[posinds]) - sum(tes[neginds])  
end

function update_ηs(test::VisitationFrequencyTest)

    VisitationFrequencyTest(
        k = test.k, l = test.l, m = test.m, n = test.n, b = test.b, 
        τ = test.τ,
        binning_summary_statistic = test.binning_summary_statistic,
        binning = test.binning,
        ηs = verified_prediction_lags(test.ηs))
end

function update_ηs(test::TransferOperatorGridTest)
    
    TransferOperatorGridTest(
        k = test.k, l = test.l, m = test.m, n = test.n, b = test.b, 
        τ = test.τ,
        binning_summary_statistic = test.binning_summary_statistic,
        binning = test.binning,
        ηs = verified_prediction_lags(test.ηs))
end

return_predictive_asymmetry(ηs, As) = As
return_predictive_asymmetry(η::Int, As) = As[1]

function predictive_asymmetry(source, target, 
        p::PredictiveAsymmetryTest{T}) where {T <: TransferEntropyCausalityTest}

    # Update the test parameters so that we have symmetric prediction lags
    test = update_ηs(p.predictive_test)
    
    # Predictions from source to target. 
    preds = causality(source, target, test)

    # Compute predictive asymmetries. The number of predictive asymmetries is half the 
    # number of prediction lags)
    ηs = test.ηs
    As = zeros(Float64, round(Int, length(ηs)/2))

    for (i, η) in enumerate(ηs[ηs .> 0])
        lag_idxs = findfirst(ηs .== -η):findfirst(ηs .== η)
        As[i] = causalbalance(ηs[lag_idxs], preds[lag_idxs])
    end

    return return_predictive_asymmetry(p.predictive_test.ηs, As)
end

function causality(source::AbstractVector{T}, target::AbstractVector{T}, p::PredictiveAsymmetryTest{CT}) where {T<:Real, CT}
    predictive_asymmetry(source, target, p)
end

# there is no need to define custom causality method for a uncertain data 
# as long as the predictive tests that it uses supports them

export 
    PredictiveAsymmetryTest, 
    predictive_asymmetry, 
    update_ηs, 
    verified_prediction_lags