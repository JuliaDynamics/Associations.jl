import DelayEmbeddings: embed
using Statistics: cor

export PairwiseAsymmetricInference, PAI
"""
    PairwiseAsymmetricInference <: CrossmapMeasure
    PairwiseAsymmetricInference(; d::Int = 2, τ::Int = -1, w::Int = 0,
        f = Statistics.cor, embed_warn = true)

The pairwise asymmetric inference (PAI) measure (McCracken & Weigel (2014)[^McCracken2014])
is a version of [`ConvergentCrossMapping`](@ref) that searches for neighbors in
*mixed* embeddings (i.e. both source and target variables included); otherwise, the
algorithms are identical.

Specifies embedding dimension `d`, embedding lag `τ` to be used, as described below,
with [`predict`](@ref) or [`crossmap`](@ref). The Theiler window `w` controls how many
temporal neighbors are excluded during neighbor searches (`w = 0` means that only the
point itself is excluded).
`f` is a function that computes the agreement between observations and
predictions (the default, `f = Statistics.cor`, gives the Pearson correlation
coefficient).

## Embedding

There are many possible ways of defining the embedding for PAI. Currently, we only
implement the *"add one non-lagged source timeseries to an embedding of the target"*
approach, which is used as an example in McCracken & Weigel's paper. Specifically:
Let `S(i)` be the source time series variable and `T(i)` be the target time series variable.
`PairwiseAsymmetricInference` produces regular embeddings with fixed dimension `d` and
embedding lag `τ` as follows:

```math
(S(i), T(i+(d-1)\\tau, \\ldots, T(i+2\\tau), T(i+\\tau), T(i)))_{i=1}^{N-(d-1)\\tau}.
```

In this joint embedding, neighbor searches are performed in the subspace spanned by
the first `D` variables, while the last variable is to be predicted.

With this convention, `τ < 0` implies "past/present values of source used to predict
target", and `τ > 0` implies "future/present values of source used to predict target".
The latter case may not be meaningful for many applications, so by default, a warning
will be given if `τ > 0` (`embed_warn = false` turns off warnings).

[^McCracken2014]:
    McCracken, J. M., & Weigel, R. S. (2014). Convergent cross-mapping and pairwise
    asymmetric inference. Physical Review E, 90(6), 062903.
"""
Base.@kwdef struct PairwiseAsymmetricInference <: CrossmapMeasure
    d::Int = 2
    τ::Int = -1
    w::Int = 0
    f::Function = cor
    embed_warn::Bool = true
end
const PAI = PairwiseAsymmetricInference

n_neighbors_simplex(measure::PairwiseAsymmetricInference) =
    (measure.d + 1) + 1 # one extra coordinate included, due to the inclusion of the target.
max_segmentlength(measure::PairwiseAsymmetricInference, x::AbstractVector) =
    length(x) - measure.d + 1
# TODO: version that takes into consideration prediction lag

function embed(measure::PairwiseAsymmetricInference, t::AbstractVector, s::AbstractVector)
    (; d, τ, w) = measure
    @assert τ != 0
    if τ > 0 && measure.embed_warn
        @warn """τ > 0. You're using future values of source to predict the target. Turn \
        off this warning by setting `embed_warn = false` in the \
        `PairwiseAsymmetricInference` constructor."""
    end
    # Convention:
    # - Negative τ := embedding vectors (s(i), t(i), t(i-1), ...), "past predicts present"
    # - Positive τ := embedding vectors (s(i), t(i), t(i+1), ...), "future predicts present"
    τs = [0; reverse(range(start=0, step=τ, stop=(d-1)*τ))]
    js = [2; repeat([1], d)]
    idxs_S̄ = 1:measure.d
    idx_t̄ = measure.d + 1 # column index of time series to be predict
    return genembed(StateSpaceSet(t, s), τs, js), idx_t̄, idxs_S̄
end
