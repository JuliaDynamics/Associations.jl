import DelayEmbeddings: embed
using Statistics: cor

export PredictiveDistanceCorrelation

"""
    PredictiveDistanceCorrelation <: CrossmapMeasure
    PredictiveDistanceCorrelation(d = 2, τ = 1, w = 0)

The predictive distance correlation (PDC) measure (Haaga, in prep),
which is used with [`predict`](@ref) and [`crossmap`](@ref).

## Arguments

- `d::Int`. The embedding dimension.
- `τ::Int`. The embedding lag.
- `w::Int`. The Theiler window, which controls how many temporal neighbors are excluded
    during neighbor searches. The default `w = 0`, which means that only the point itself
    is excluded.

## Description

`PredictiveDistanceCorrelation` is a multivariate extension to [`ConvergentCrossMapping`](@ref) and
[`PairwiseAsymmetricInference`](@ref), where agreement between predictions and
and observed values are measured using [`distance_correlation`](@ref).

Assume we want to predict a `target` StateSpaceSet using
nearest-neighbor information from a `source` StateSpaceSet. The PDC algorithm proceeds as
follows for each `i ∈ 1, 2, …, N`, where
`N = length(target) = length(source) = length(cond)`.

    1. Let `T̄ₛ` be an `N`-element vector, where we will store predictions.
    2. Locate `sᵢ ∈ source` and find its `Dt + 1` nearest neighbors. Label the
    indices `n1, n2, …, nK`, where `K = Dt + 1`.
    3. Compute the prediction `t̄ᵢ ∈ ℝᴰᵗ` using some linear combination of the points
        `target[n1], target[n2], …, target[nK]` (which form an enclosing simplex
        around the point `tᵢ ∈ target`). Weights for the linear combination are
        determined by the measure `m`.
    4. Set `T̄ₛ[i] = t̄ᵢ`.

Now, both observations `tᵢ ∈ ℝᴰᵗ` and predictions `t̄ᵢ ∈ ℝᴰᵗ` for `i ∈ 1, 2, …, N`.
Finally, compute the degree of agreement between observations and predictions as the
distance correlation `ρT̄ₛT = dcov(T̄ₛ, target)`. We have `ρT̄ₛT ∈ [0, 1]`, `ρT̄ₛT = 0`
is no agreement and `ρT̄ₛT = 1` is perfect agreement.
"""
Base.@kwdef struct PredictiveDistanceCorrelation <: CrossmapMeasure
    d::Int = 2
    τ::Int = 1
    w::Int = 0 # Theiler window
end

function crossmap(m::PredictiveDistanceCorrelation, T::AbstractStateSpaceSet, S::AbstractStateSpaceSet, C::AbstractStateSpaceSet)
    TC = StateSpaceSet(T, C)
    T̄C̄ₛ = predict(m, StateSpaceSet(TC), S);
    ρT̄C̄ₛTC = distance_correlation(T̄C̄ₛ, TC)
    return ρT̄C̄ₛTC
end

function crossmap(m::PredictiveDistanceCorrelation, T::AbstractStateSpaceSet, S::AbstractStateSpaceSet)
    T̄ₛ = predict(m, T, S);
    ρT̄ₛT = distance_correlation(T̄ₛ, T)
    return ρT̄ₛT
end

function predict(measure::PredictiveDistanceCorrelation, T̄::AbstractStateSpaceSet{DT}, S̄::AbstractStateSpaceSet{DS}) where {DT, DS}
    @assert length(S̄) == length(T̄)
    (; d, τ, w) = measure
    # Tree structure must be created for every L, because we can't include data
    # outside the considered time range.

    # The number of neighbors depend on the type of cross map measure. We could make
    # this a tunable parameter, but for now, just stick with dim(embedding) + 1.
    nnd = dimension(T̄) + 1
    tree = KDTree(S̄, Euclidean())
    # Todo: maybe re-use the same tree, but with a more elaborate skip function?
    # Not sure what is fastest. Need to experiment...
    nnidxs, ds = bulksearch(tree, S̄, NeighborNumber(nnd), Theiler(w))
    T̂ = Vector{SVector{DT}}(undef, 0) # predicted values
    u = zeros(MVector{nnd}) # one extra neighbor due to the extra coordinate from t
    w = zeros(MVector{nnd}) # one extra neighbor due to the extra coordinate from t

    t̂ = zeros(MVector{DT})# prediction vector.
    for (i, (nnidxsᵢ, dᵢ)) in enumerate(zip(nnidxs, ds))
        u .= exp.(-dᵢ ./ dᵢ[1])
        w .= u ./ sum(u)
        # The predicted vector t̂ is a linear combination of the points in T̄[nnidxsᵢ], where
        # weights are determines by neighbor distances to point sᵢ ∈ S̄
        t̂ .= 0 # re-zero
        for d = 1:DT
            t̂ .+= w[d] .* T̄[nnidxsᵢ][d]
        end
        push!(T̂, t̂)
    end
    return StateSpaceSet(T̂)
end


export condmap
# Experimental: conditional prediction M(S → T | C) (Haaga et al, paper in prep)
function condmap(m::CrossmapMeasure, t::AbstractStateSpaceSet, s::AbstractStateSpaceSet,
        c::AbstractStateSpaceSet)
    T = StateSpaceSet(t)
    S = StateSpaceSet(s)
    C = StateSpaceSet(c)
    TC = StateSpaceSet(T, C)

    # Unconditional and conditional cross map, and their distance
    # correlations.
    T̄ₛ = predict(m, T, S);
    T̄C̄ₛ = predict(m, StateSpaceSet(TC), S);
    ρT̄ₛT = distance_correlation(T̄ₛ, T)
    ρT̄C̄ₛTC = distance_correlation(T̄C̄ₛ, TC)

    # If predictions T̄C̄ₛ are better than T̄ₛ, then we take that as evidence that
    #
    # If including information about C in the target increases
    # quality of predictions, then C also contributes to S, i.e.
    # both X and Z contribute to Y.
    ΔρTS_given_C = ρT̄C̄ₛTC - ρT̄ₛT
    return ΔρTS_given_C
end



# ## Arguments

# - *`m::CrossmapMeasure`*. A cross-map measure, e.g. [`ConvergentCrossMapping`](@ref), that determines
#     the procedure for cross mapping. *Note* `crossmap` doesn't embed input data,
#     so the embedding parameters of `m` are ignored.
# - *`target::AbstractVector`*. A scalar-valued time series.
# - *`source::AbstractStateSpaceSet`*. A `Ds`-dimensional source StateSpaceSet.

# ## Description

# The following procedure is repeated for each `i ∈ 1, 2, …, N`, where
# `N = length(target) = length(source) = length(cond)`. Let `T̄ₛ` be an
# `N`-element vector, where we will store predictions.

# 1. Locate `sᵢ ∈ source` and find its `Dt + 1` nearest neighbors. Label the
# indices `n1, n2, …, nK`, where `K = Dt + 1`.
# 2. Compute the prediction `t̄ᵢ ∈ ℝᴰᵗ` using some linear combination of the points
#     `target[n1], target[n2], …, target[nK]` (which form an enclosing simplex
#     around the point `tᵢ ∈ target`). Weights for the linear combination are
#     determined by the measure `m`. Set `T̄ₛ[i] = t̄ᵢ`.

# Now all `tᵢ ∈ ℝᴰᵗ` and all `t̄ᵢ ∈ ℝᴰᵗ`. We then compute the degree of agreement between
# observations and predictions as the distance correlation `ρT̄ₛT = dcov(T̄ₛ, target)`,
# `ρT̄ₛT ∈ [0, 1]`, and `ρT̄ₛT = 0` is no agreement and `ρT̄ₛT = 1` is perfect agreement.
