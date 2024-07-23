import DelayEmbeddings: embed
using DelayEmbeddings: genembed
using Statistics: cor

export ConvergentCrossMapping, CCM
"""
    ConvergentCrossMapping <: CrossmapMeasure
    ConvergentCrossMapping(; d::Int = 2, τ::Int = -1, w::Int = 0,
        f = Statistics.cor, embed_warn = true)

The convergent cross mapping measure [Sugihara2012](@cite).

## Usage

- Use with [`association`](@ref) together with a [`CrossmapEstimator`](@ref) to compute the 
    cross-map correlation between input variables.

```julia
est_single_lib = RandomSegment(ConvergentCrossMapping(); libsizes = 100)
est_multiple_libs = RandomSegment(ConvergentCrossMapping; libsizes = [100, 200, 300])
association(est_single_lib, x, y) → ρ ∈ [-1, 1]
association(est_multiple_libs, x, y) → Vector{ρᵢ} (ρᵢ ∈ [-1, 1])
```

## Compatible estimators

- [`RandomSegment`](@ref)
- [`RandomVectors`](@ref)
- [`ExpandingSegment`](@ref)

## Description

Specifies embedding dimension `d`, embedding lag `τ` to be used, as described below,
with [`predict`](@ref) or [`crossmap`](@ref). The Theiler window `w` controls how many
temporal neighbors are excluded during neighbor searches (`w = 0` means that only the
point itself is excluded).
`f` is a function that computes the agreement between observations and
predictions (the default, `f = Statistics.cor`, gives the Pearson correlation
coefficient).

## Embedding

Let `S(i)` be the source time series variable and `T(i)` be the target time series variable.
This version produces regular embeddings with fixed dimension `d` and embedding lag
`τ` as follows:

```math
( S(i), S(i+\\tau), S(i+2\\tau), \\ldots, S(i+(d-1)\\tau, T(i))_{i=1}^{N-(d-1)\\tau}.
```

In this joint embedding, neighbor searches are performed in the subspace spanned by
the first `D-1` variables, while the last (`D`-th) variable is to be predicted.

With this convention, `τ < 0` implies "past/present values of source used to predict
target", and `τ > 0` implies "future/present values of source used to predict target".
The latter case may not be meaningful for many applications, so by default, a warning
will be given if `τ > 0` (`embed_warn = false` turns off warnings).

## Estimation

- [Example 1](@ref example_ConvergentCrossMapping_RandomVectors). 
    Estimation with [`RandomVectors`](@ref) estimator.
- [Example 2](@ref example_ConvergentCrossMapping_RandomSegment). 
    Estimation with [`RandomSegment`](@ref) estimator.
- [Example 3](@ref example_ConvergentCrossMapping_reproducing_sugihara): Reproducing 
    figures from [Sugihara2012](@citet).
"""
Base.@kwdef struct ConvergentCrossMapping <: CrossmapMeasure
    d::Int = 2
    τ::Int = -1
    w::Int = 0
    f::Function = cor
    embed_warn::Bool = true
end
const CCM = ConvergentCrossMapping

n_neighbors_simplex(definition::ConvergentCrossMapping) = definition.d + 1
max_segmentlength(definition::ConvergentCrossMapping, x::AbstractVector) =
    length(x) - definition.d + 1
# TODO: version that takes into consideration prediction lag

function embed(definition::ConvergentCrossMapping, t::AbstractVector, s::AbstractVector)
    (; d, τ, w, f) = measure
    if τ > 0 && definition.embed_warn
        @warn """τ > 0. You're using future values of source to predict the target. Turn \
        off this warning by setting `embed_warn = false` in the \
        `PairwiseAsymmetricInference` constructor."""
    end
    @assert τ != 0
    # Convention:
    # - Negative τ := embedding vectors (s(i), s(i-1), ..., t(i)), "past predicts present"
    # - Positive τ := embedding vectors (s(i), s(i+1), ..., t(i)), "future predicts present"
    τs = [0:τ:(d-1)*τ; 0]
    js = [repeat([1], d); 2]
    idxs_S̄ = 1:length(js) - 1
    idx_t̄ = length(js)# column index of time series to be predict
    genembed(StateSpaceSet(s, t), τs, js), idx_t̄, idxs_S̄
end
