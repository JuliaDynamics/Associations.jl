import DelayEmbeddings: embed
using Statistics: cor

export PairwiseAsymmetricInference
"""
    PairwiseAsymmetricInference <: CrossmapMeasure
    PairwiseAsymmetricInference(d = 2, τ = 1, w = 0, f = Statistics.cor)

The pairwise asymmetric inference (PAI) measure, as given in McCracken &
Weigel (2014)[^McCracken2014].

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
`PairwiseAsymmetricInference` produces regular embeddings with fixed dimension `d` and embedding lag
`τ` as follows:

```math
( S(i), X(i-\\tau), S(i-2\\tau), \\ldots, S(i-(d-1)\\tau, T(i)))_{i=1}^{N-(d-1)\\tau}.
```

Note the difference from [`ConvergentCrossMapping`](@ref): for `PairwiseAsymmetricInference`, we include an additional
non-lagged version of the *target* time series in the embedding.

[^McCracken2014]:
    McCracken, J. M., & Weigel, R. S. (2014). Convergent cross-mapping and pairwise
    asymmetric inference. Physical Review E, 90(6), 062903.
"""
Base.@kwdef struct PairwiseAsymmetricInference <: CrossmapMeasure
    d::Int = 2
    τ::Int = 1
    w::Int = 0
    f::Function = Statistics.cor
end

n_neighbors_simplex(measure::PairwiseAsymmetricInference) = (measure.d + 1) + 1 # one extra coordinate included
max_segmentlength(x::AbstractVector, measure::PairwiseAsymmetricInference) = length(x) - measure.d + 1

function embed(measure::PairwiseAsymmetricInference, s::AbstractVector, t::AbstractVector)
    (; d, τ, w) = measure
    τs = [0:τ:(d-1)*τ; 0; 0]
    js = [repeat([1], d); 2; 2]
    genembed(Dataset(s, t), τs, js)
end
