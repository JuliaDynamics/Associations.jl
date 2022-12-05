import DelayEmbeddings: embed
using Statistics: cor

export ConvergentCrossMapping
"""
    ConvergentCrossMapping <: CrossmapMeasure
    ConvergentCrossMapping(d = 2, τ = 1, w = 0, f = Statistics.cor)

The convergent cross mapping (CCM) measure, as given in Sugihara et al.
(2012)[^Sugihara2012].

## Arguments

- `d::Int`. The embedding dimension.
- `τ::Int`. The embedding lag.
- `f::Function`. A function that computes the agreement between observations and
    predictions. The default, `f = Statistics.cor`, gives the Pearson correlation
    coefficient.
- `w::Int`. The Theiler window, which controls how many temporal neighbors are excluded
    during neighbor searches. The default `w = 0`, which means that only the point itself
    is excluded.

## Description

Let `S(i)` be the source time series variable and `T(i)` be the target time series variable.
This version produces regular embeddings with fixed dimension `d` and embedding lag
`τ` as follows:

```math
( S(i), S(i-\\tau), S(i-2\\tau), \\ldots, S(i-(d-1)\\tau))_{i=1}^{N-(d-1)\\tau}.
```

[^Sugihara2012]:
    Sugihara, G., May, R., Ye, H., Hsieh, C. H., Deyle, E., Fogarty, M., & Munch, S.
    (2012). Detecting causality in complex ecosystems. science, 338(6106), 496-500.
"""
Base.@kwdef struct ConvergentCrossMapping <: CrossmapMeasure
    d::Int = 2
    τ::Int = 1
    w::Int = 0
    f::Function = Statistics.cor
end

n_neighbors_simplex(measure::ConvergentCrossMapping) = measure.d + 1
max_segmentlength(x::AbstractVector, measure::ConvergentCrossMapping) = length(x) - measure.d + 1

function embed_for_crossmap(measure::ConvergentCrossMapping, t::AbstractVector, s::AbstractVector)
    (; d, τ, w) = measure
    τs = [0:τ:(d-1)*τ; 0]
    js = [repeat([1], d); 2]
    genembed(Dataset(s, t), τs, js)
end

function estimate(measure::ConvergentCrossMapping, est::ExpandingSegment, t::AbstractVector, s::AbstractVector)
    (; d, τ, w) = measure
    (length(x) - d) > 10 + d + 1 || error("Input time series are too short")

    # Embedding both time series together ensures that we don't mess up indices later.
    # The first `d` variables of the embedding are for x, the last is for y
    τs = [0:τ:(d-1)*τ; 0]
    js = [repeat([1], d); 2]

    # Embed together, so we don't mess up time indices
    mixed_embedding = genembed(Dataset(x, y), τs, js)
    S̄ = mixed_embedding[:, 1:d]
    target = mixed_embedding[:, end]
    L = max_segmentlength(x, measure)
    N = length(S̄)
    Ls = max((L ÷ 10), 10):10:L

    ρs = zeros(length(Ls))
    for (i, L) in enumerate(Ls)
        s = S̄[1:L]
        t = target[1:L]
        t̂ = predict(measure, s, t)
        ρ[i] = cor(t, t̂)
    end
    return Ls, ρs
end
