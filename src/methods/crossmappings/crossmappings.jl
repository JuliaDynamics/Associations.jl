include("utils.jl")
import DelayEmbeddings: embed
using Neighborhood: Euclidean, KDTree, NeighborNumber, Theiler
using Neighborhood: bulksearch
using StaticArrays: MVector

export predict, crossmap
export CrossmapMeasure, CrossmapEstimator
export ExpandingSegment, RandomSegment

"""
The supertype for all cross-map measures.
"""
abstract type CrossmapMeasure end

"""
The supertype for all cross-map estimators.
"""
abstract type CrossmapEstimator end

segment_length_error() = "Segment lengths can be inferred only if both a cross-map " *
    "measure and an input time series is provided. " *
    "Do e.g. `ExpandingSegment(CCM(), x)`, where `x` is some time series."

"""
    max_segmentlength(x::AbstractVector, measure::CrossmapMeasure)

Given an input vector `x`, which is representative of the size of the other input vectors
too, compute the maximum segment/library length that can be used for predictions.
"""
function max_segmentlength end

"""
    embed(measure::CrossmapMeasure, target::AbstractVector,
        sources::AbstractVector...) → emb::Dataset{D}

Jointly embed the input vector `target`, together with one or more vectors `s ∈ sources`,
according to the given [`CrossmapMeasure`](@ref).

This produces `emb`, a `D`-dimensional `Dataset` where

- The last column is always the non-lagged `target` variable. Typically, this is
    the variable we're trying to predict.
- The `D-1` first columns are the (non)lagged versions of each source time series
    `s ∈ sources`. Typically, `emb[:, 1:D-1]` is the subspace in which neighborhood
    searches are done, which forms the basis of cross-map [`predict`](@ref)ions.
"""
function embed end

"""
    crossmap(measure, target, source::AbstractDataset) → ρt̄ₛt

Like [`predict`](@ref), but return the correspondence between `target` and `t̄ₛ`,
according to the given cross-map `measure`.

See also: [`CrossmapMeasure`](@ref), [`ConvergentCrossMapping`](@ref), [`PairwiseAsymmetricInference`](@ref).
"""
function crossmap end

function crossmap(measure::CrossmapMeasure, target::AbstractVector, source::AbstractDataset)
    t̄ₛ = predict(measure, target, source);
    ρt̄ₛt = measure.f(t̄ₛ, target)
    return ρt̄ₛt
end

"""
    predict(measure, target, source::AbstractDataset{Ds}) → t̂

For each `i ∈ 1, 2, …, N = length(target) = length(source)`, produce the cross-map
prediction `t̂ₛ[i]` by a linear combination of points in `target`, where the selection of
points and weights are determined by the `D+1` nearest neighbors of the point `source[i]`.

## Returns

Returns a vector of predictions `t̂`, where `t̂ₛ[i]` is the prediction for `target[i]`.

You can either compute the agreement between predictions and observations manually,
or use [`crossmap`](@ref) instead, which directly returns the correspondence between
`target` and `t̂ₛ` according to `measure`.

## Arguments

- *`measure::CrossmapMeasure`*. A cross-map measure, e.g. [`ConvergentCrossMapping`](@ref), which
    determines which cross-map procedure is applied. *Note* `crossmap` doesn't embed
    input data, so any embedding parameters in `measure` are ignored.
- *`target::AbstractVector`*. A scalar-valued vector.
- *`source::AbstractDataset`*. A `Ds`-dimensional source dataset.

## Data requirements

Assumes that `target` and `source` have already been jointly embedded using [`embed`](@ref),
so that indexing is correct.

See also: [`CrossmapMeasure`](@ref), [`ConvergentCrossMapping`](@ref), [`PairwiseAsymmetricInference`](@ref).
"""
function predict end

# """
#     predict(measure::CrossmapMeasure, t::AbstractVector, S̄::AbstractDataset) → t̂
#     predict(measure::CrossmapMeasure, T̄::AbstractDataset, S̄::AbstractDataset) → t̂

# A version of `predict` that uses the already-constructed embedding `S̄` to predict
# univariate timeseries `t` using the given `measure`.

# To obtain the correlation between observed and predicted values, you can do:
# - `cor(t, t̂)` if `t isa AbstractVector`
# - `fastcor(t, t̂)` if `t isa Dataset`. `fastcor` uses canonical correlation analysis
#      to project t̂ and S̄ into a common 2D space where `cor(t, t̂)` to compute the
#      correlation between observed and predicted values in the reduced 2D space.
# """
function predict(measure::CrossmapMeasure, t::AbstractVector, S̄::AbstractDataset)
    @assert length(S̄) == length(t)
    (; d, τ, w) = measure
    # Tree structure must be created for every L, because we can't include data
    # outside the considered time range.

    # The number of neighbors depend on the type of cross map measure. We could make
    # this a tunable parameter, but for now, just stick with dim(embedding) + 1.
    nnd = dimension(S̄) + 1
    tree = KDTree(S̄, Euclidean())
    # Todo: maybe re-use the same tree, but with a more elaborate skip function?
    # Not sure what is fastest. Need to experiment...
    nnidxs, ds = bulksearch(tree, S̄, NeighborNumber(nnd), Theiler(w))
    t̂ = zeros(length(S̄))
    u = zeros(MVector{nnd}) # one extra neighbor due to the extra coordinate from t
    w = zeros(MVector{nnd}) # one extra neighbor due to the extra coordinate from t
    for (i, (nnidxsᵢ, dᵢ)) in enumerate(zip(nnidxs, ds))
        u .= exp.(-dᵢ ./ dᵢ[1])
        w .= u ./ sum(u)
        # Predict using weights computed from source `s` applied to values from target `t`
        t̂[i] = sum(w .* t[nnidxsᵢ])
    end
    return t̂
end

# Multivariate version. Used by PredictiveDistanceCorrelation.
function predict(measure::CrossmapMeasure, T̄::AbstractDataset{DT},
        S̄::AbstractDataset{DS}) where {DT, DS}
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
        # The predicted vector t̂ is the center of mass, weighted by distances from S̄
        t̂ .= 0 # re-zero
        for d = 1:DT
            t̂ .+= w[d] .* T̄[nnidxsᵢ][d]
        end
        t̂ ./= nnd
        push!(T̂, t̂)
    end
    return Dataset(T̂)
end

include("estimators/estimators.jl")
include("measures/measures.jl")
