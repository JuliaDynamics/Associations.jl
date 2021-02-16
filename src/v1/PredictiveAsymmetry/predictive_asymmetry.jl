
import TransferEntropy: TransferEntropyEstimator, EmbeddingTE, transferentropy
export predictive_asymmetry

function causalbalance(lags, tes)
    posinds = findall(lags .> 0)
    neginds = findall(lags .< 0)
    sum(tes[posinds]) - sum(tes[neginds])  
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
    predictive_asymmetry(source, target, [cond], 
        estimator::TransferEntropyEstimator, ηs; 
        d𝒯 = 1, dT = 1, dS = 1, τT = -1, τS = -1, 
        [dC = 1, τC = -1,],
        normalize::Bool = false, f::Real = 1.0)

Compute the predictive asymmetry[^Haaga2020] 𝔸(`source` → `target`) over prediction lags 
`ηs`, using the given transfer entropy `estimator` and embedding parameters `d𝒯`, `dT`, 
`dS`, `τT`, `τS`.

If `cond` is provided, compute 𝔸(`source` → `target` | `cond`). Then, `dC` and `τC` controls 
the embedding dimension and embedding lag for the conditional variable.

## Normalization (hypothesis test)

If `normalize == true` (the default), then compute the normalized predictive asymmetry 𝒜. 

In this case, for each ``\\eta`` in `ηs`, compute 𝒜(η) by normalizing 𝔸(η) to some fraction `f` of the 
mean transfer entropy over prediction lags ``-\\eta, ..., \\eta`` (exluding lag 0). 

Haaga et al. (2020)[^Haaga2020] uses a normalization with `f=1.0` as a built-in hypothesis test, 
avoiding more computationally costly surrogate testing.

## Examples

It is recommended to use either the rectangular binning-based methods or the symbolic estimators 
for the fastest computations. 

```julia
# Some example time series
x, y, z = rand(100), rand(100), rand(100)

# Define prediction lags and estimation method
ηs = 1:5
method = VisitationFrequency(RectangularBinning(5))

# 𝔸(x → y) and  𝔸(x → y | z)
𝔸reg  = predictive_asymmetry(x, y, method, ηs, normalize = false)
𝔸cond = predictive_asymmetry(x, y, z, method, ηs, normalize = false)

# 𝒜(x → y) and 𝒜(x → y | z), using different normalization factors
𝒜reg = predictive_asymmetry(x, y, ηs, method, f = 1.0) # normalize == true by default
𝒜cond = predictive_asymmetry(x, y, z, ηs, method, f = 1.5) # normalize == true by default
```

For the symbolic estimators, make sure that the maximum prediction lag η stays 
small. This is because the symbolization procedure uses delay embedding vectors 
of dimension `m` if the motif length is `m` (so the actual maximum prediction lag 
used is `maximum(ηs)*m`, which may be too large if `maximum(ηs)` is too large).


```julia
# Some example time series
x, y, z = rand(100), rand(100), rand(100)

# Define prediction lags and estimation method
ηs = 1:3 # small prediction lags
estimator = VisitationFrequency(RectangularBinning(4))

# 𝒜(x → y)
predictive_asymmetry(x, y, estimator, ηs, normalize = true) 

# 𝒜(x → y | z)
predictive_asymmetry(x, y, z, estimator, ηs, normalize = true) 
```

[^Haaga2020]: Haaga, Kristian Agasøster, David Diego, Jo Brendryen, and Bjarte Hannisdal. "A simple test for causality in complex systems." arXiv preprint arXiv:2005.01860 (2020).
"""
function predictive_asymmetry end

function predictive_asymmetry(source, target, estimator, ηs; 
        normalize = false, f::Real = 1.0,
        d𝒯 = 1, dT = 1, dS = 1, τT = -1, τS = -1)
    
    Nη = length(ηs)
    te_fws = zeros(Nη)
    te_bws = zeros(Nη)
    𝔸s = zeros(Nη)

    for (i, η) in enumerate(ηs)
        te_fws[i] = transferentropy(source, target, estimator, d𝒯 = d𝒯, dT = dT, dS = dS, τT = τT, τS = τS, η𝒯 = η)
        te_bws[i] = transferentropy(source, target, estimator, d𝒯 = d𝒯, dT = dT, dS = dS, τT = τT, τS = τS, η𝒯 = -η)
        
        𝔸s[i] = sum(te_fws[1:i]) - sum(te_bws[1:i])
    end
    
    if normalize
        for (i, η) in enumerate(ηs)
            𝔸ᵢ = sum(te_fws[1:i]) - sum(te_bws[1:i])
            avg_te = 1/(2*η) * (sum(te_fws[1:i]) + sum(te_bws[1:i]))
            𝔸s[i] = 𝔸ᵢ/(f*avg_te)
        end
    end
    
    return 𝔸s
end

function predictive_asymmetry(source, target, cond, estimator, ηs; 
        normalize = true, f::Real = 1.0,
        d𝒯 = 1, dT = 1, dS = 1, dC = 1, τT = -1, τS = -1, τC = -1)
    
    Nη = length(ηs)
    all(ηs .> 0) || throw(ArgumentError("all ηs must be >= 1, got $(ηs)"))
    issorted(ηs) || throw(ArgumentError("ηs must be provided in increasing order, got $(ηs)"))

    te_fws = zeros(Nη)
    te_bws = zeros(Nη)
    𝔸s = zeros(Nη)
    
    for (i, η) in enumerate(ηs)
        te_fws[i] = transferentropy(source, target, cond, estimator, d𝒯 = d𝒯, dT = dT, dS = dS, τT = τT, τS = τS, dC = dC, τC = τC, η𝒯 = η)
        te_bws[i] = transferentropy(source, target, cond, estimator, d𝒯 = d𝒯, dT = dT, dS = dS, τT = τT, τS = τS, dC = dC, τC = τC, η𝒯 = -η)
        
        𝔸s[i] = sum(te_fws[1:i]) - sum(te_bws[1:i])

        if normalize 
            avg_te = 1/(2*η) * (sum(te_fws[1:i]) + sum(te_bws[1:i]))
            𝔸s[i] = 𝔸s[i]/(f*avg_te) # after normalization, this is now 𝒜s(i)
        end
    end
    
    return 𝔸s
end