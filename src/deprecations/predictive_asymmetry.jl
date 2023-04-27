export predictive_asymmetry

function causalbalance(lags, tes)
    posinds = findall(lags .> 0)
    neginds = findall(lags .< 0)
    sum(tes[posinds]) - sum(tes[neginds])
end

function verified_prediction_lags(lag::Int)
    # If only one η is provided (either positive or negative), just make a prediction
    # for that lag.
    ηs = [-abs(lag), abs(lag)]
end

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

"""
    predictive_asymmetry(estimator::TransferEntropyEstimator, ηs; s, t, [c],
        dTf = 1, dT = 1, dS = 1, τT = -1, τS = -1, [dC = 1, τC = -1],
        normalize::Bool = false, f::Real = 1.0, base = 2) → Vector{Float64}

Compute the predictive asymmetry[^Haaga2020] 𝔸(`s` → `t`) for source time series `s` and
target time series `t` over prediction lags `ηs`, using the given `estimator` and embedding
parameters `dTf`, `dT`, `dS`, `τT`, `τS` (see also [`EmbeddingTE`](@ref))

If a conditional time series `c` is provided, compute 𝔸(`s` → `t` | `c`). Then, `dC` and
`τC` controls the embedding dimension and embedding lag for the conditional variable.

## Returns

Returns a vector containing the predictive asymmetry for each value of `ηs`.

## Normalization (hypothesis test)

If `normalize == true` (the default), then compute the normalized predictive asymmetry 𝒜.
In this case, for each ``\\eta`` in `ηs`, compute 𝒜(η) by normalizing 𝔸(η) to some fraction `f` of the
mean transfer entropy over prediction lags ``-\\eta, ..., \\eta`` (exluding lag 0).
Haaga et al. (2020)[^Haaga2020] uses a normalization with `f=1.0` as a built-in hypothesis test,
avoiding more computationally costly surrogate testing.

## Estimators

Any [estimator](@ref) that works for [`transferentropy`](@ref) will also work with
`predictive_asymmetry`. Check the online documentation for compatiable estimators.

## Examples

```julia
using CausalityTools
# Some example time series
x, y = rand(100), rand(100)
# 𝔸(x → y) over prediction lags 1:5
𝔸reg  = predictive_asymmetry(x, y, VisitationFrequency(RectangularBinning(5)), 1:5)
```

!!! info "Experimental!"
    This is a method that does not yet appear in a peer-reviewed scientific journal.
    Feel free to use, but consider it experimental for now. It will reappear in
    a 2.X release in new form once published in a peer-reviewed journal.

[^Haaga2020]:
    Haaga, Kristian Agasøster, David Diego, Jo Brendryen, and Bjarte Hannisdal.
    "A simple test for causality in complex systems."
    arXiv preprint arXiv:2005.01860 (2020).
"""
function predictive_asymmetry end

function check_ηs(ηs)
    all(ηs .> 0) || throw(ArgumentError("all ηs must be >= 1, got $(ηs)"))
    issorted(ηs) || throw(ArgumentError("ηs must be provided in increasing order, got $(ηs)"))
end

const VALID_PA_ESTIMATORS = Union{
    ProbabilitiesEstimator,
    DifferentialEntropyEstimator,
    MutualInformationEstimator,
    ConditionalMutualInformationEstimator,
    TransferEntropyEstimator
    }

Base.@deprecate predictive_asymmetry(source::AbstractVector, target::AbstractVector, estimator::VALID_PA_ESTIMATORS, ηs;
    normalize = false, f = 1.0, d𝒯 = 1, dT = 1, dS = 1, τT = -1, τS = -1) predictive_asymmetry(estimator::VALID_PA_ESTIMATORS, ηs, source, target;
        normalize = false, f = 1.0, dTf = 1, dT = 1, dS = 1, τT = -1, τS = -1, base = 2) false

Base.@deprecate predictive_asymmetry(source, target, cond, estimator::VALID_PA_ESTIMATORS, ηs;
    normalize = false, f = 1.0,
    d𝒯 = 1, dT = 1, dS = 1, dC = 1, τT = -1, τS = -1, τC = -1) predictive_asymmetry(estimator::VALID_PA_ESTIMATORS, ηs, source, target, cond;
    normalize = false, f = 1.0,
    dTf = 1, dT = 1, dS = 1, dC = 1, τT = -1, τS = -1, τC = -1, base = 2) false

function predictive_asymmetry(estimator::VALID_PA_ESTIMATORS, ηs, source, target;
        normalize = false, f::Real = 1.0,
        dTf = 1, dT = 1, dS = 1, τT = -1, τS = -1, base = 2)
    check_ηs(ηs)
    Nη = length(ηs)

    te_fws = zeros(Nη)
    te_bws = zeros(Nη)
    𝔸s = zeros(Nη)

    for (i, η) in enumerate(ηs)
        emb_fw = EmbeddingTE(dTf = dTf, dT = dT, dS = dS, τT = τT, τS = τS, ηTf = η)
        emb_bw = EmbeddingTE(dTf = dTf, dT = dT, dS = dS, τT = τT, τS = τS, ηTf = -η)
        te_fw = TEShannon(; base, embedding = emb_fw)
        te_bw = TEShannon(; base, embedding = emb_bw)
        te_fws[i] = transferentropy(te_fw, estimator, source, target)
        te_bws[i] = transferentropy(te_bw, estimator, source, target)

        if normalize
            𝔸ᵢ = (sum(te_fws[1:i]) - sum(te_bws[1:i])) / η
            avg_te = (sum(te_fws[1:i]) + sum(te_bws[1:i])) / (2*η)
            𝔸s[i] = 𝔸ᵢ / (f*avg_te)
        else
            𝔸ᵢ = (sum(te_fws[1:i]) - sum(te_bws[1:i])) / η
            𝔸s[i] = 𝔸ᵢ
        end
    end

    return 𝔸s
end

function predictive_asymmetry(estimator::VALID_PA_ESTIMATORS, ηs, x::AbstractVector{<:Real}...;
        normalize = false, f::Real = 1.0,
        dTf = 1, dT = 1, dS = 1, dC = 1, τT = -1, τS = -1, τC = -1, base = 2)
    @warn """`predictive_asymmetry` will be deprecated in CausalityTools 2.0, \
    and moved to the "Experimental" section until it appears in a peer-reviewed journal. \
    It will reappear with a new signature in a 2.X release."""
    check_ηs(ηs)
    Nη = length(ηs)

    te_fws = zeros(Nη)
    te_bws = zeros(Nη)
    𝔸s = zeros(Nη)

    for (i, η) in enumerate(ηs)
        emb_fw = EmbeddingTE(dTf = dTf, dT = dT, dS = dS, τT = τT, τS = τS, ηTf = η)
        emb_bw = EmbeddingTE(dTf = dTf, dT = dT, dS = dS, τT = τT, τS = τS, ηTf = -η)
        te_fw = TEShannon(; base, embedding = emb_fw)
        te_bw = TEShannon(; base, embedding = emb_bw)
        te_fws[i] = transferentropy(te_fw, estimator, x...)
        te_bws[i] = transferentropy(te_bw, estimator, x...)
        if normalize
            𝔸ᵢ = (sum(te_fws[1:i]) - sum(te_bws[1:i])) / η
            avg_te = (sum(te_fws[1:i]) + sum(te_bws[1:i])) / (2*η)
            𝔸s[i] = 𝔸ᵢ / (f*avg_te)
        else
            𝔸ᵢ = (sum(te_fws[1:i]) - sum(te_bws[1:i])) / η
            𝔸s[i] = 𝔸ᵢ
        end
    end

    return 𝔸s
end
