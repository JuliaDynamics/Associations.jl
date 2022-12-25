include("utils.jl")
import DelayEmbeddings: embed
export embed
using Neighborhood: Euclidean, KDTree, NeighborNumber, Theiler
using Neighborhood: bulksearch
using StaticArrays: MVector
using StateSpaceSets: AbstractDataset

export predict
export crossmap
export CrossmapMeasure
export CrossmapEstimator
export CrossmapEnsemble

"""
The supertype for all cross-map measures

Currently implemented measures are:
- [`ConvergentCrossMapping`](@ref), or [`CCM`](@ref) for short.
- [`PairwiseAsymmetricInference`](@ref), or [`PAI`](@ref) for short.
"""
abstract type CrossmapMeasure end

"""
    CrossmapEstimator{LIBSIZES, RNG}

A parametric supertype for all cross-map estimators, which are used with [`predict`](@ref) and
[`crossmap`](@ref).

Because the type of the library may differ between estimators, and because RNGs from
different packages may be used, subtypes must implement the `LIBSIZES` and `RNG`
type parameters.

For efficiency purposes, subtypes may contain mutable containers that can be re-used
for ensemble analysis (see [`CrossmapEnsemble`](@ref)).

## Libraries

A cross-map estimator uses the concept of "libraries". A library is essentially just
a reference to a set of points, and usually, a library refers to *indices* of points,
not the actual points themselves.

For example, for timeseries, `RandomVectors(libsizes = 50:25:100)` produces three
separate libraries, where the first contains 50 randomly selected time indices,
the second contains 75 randomly selected time indices, and the third contains
100 randomly selected time indices. This of course assumes that all quantities involved
can be indexed using the same time indices, meaning that the concept of "library"
only makes sense *after* relevant quantities have been *jointly* embedded, so that they
can be jointly indexed. For non-instantaneous prediction, the maximum possible library
size shrinks with the magnitude of the index/time-offset for the prediction.

For spatial analyses (not yet implemented), indices could be more complex and involve
multi-indices.
"""
abstract type CrossmapEstimator{LIBSIZES, RNG} end

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
function embed(measure::CrossmapMeasure, args...) end

"""
    crossmap(measure::CrossmapMeasure, t̄::AbstractVector, S̄::AbstractDataset) → ρ
    crossmap(measure::CrossmapMeasure, target::AbstractVector, source::AbstractVector) → ρ

Compute the cross map estimates between time-aligned time series `t̄` and
source embedding `S̄`, or between raw time series `t` and `s`.

This is just a wrapper around [`predict`](@ref) that simply returns the correspondence
measure between the source and the target.
"""
function crossmap end

# Generic implementation. This is the ConvergentCrossMapping/PairwiseAsymmetricInference
# approach, but other `CrossmapMeasure`s may do something different. Custom
# implementations go in a relevant `measures/CustomMeasure.jl` file.

"""
    predict(measure::CrossmapMeasure, t̄::AbstractVector, S̄::AbstractDataset) → t̂ₛ
    predict(measure::CrossmapMeasure, target::AbstractVector, source::AbstractVector) → t̂ₛ, t̄, ρ

Perform point-wise cross mappings between source embeddings and target time series
according to the algorithm specified by the given cross-map `measure` (e.g.
[`ConvergentCrossMapping`](@ref) or [`PairwiseAsymmetricInference`](@ref)).

- **First method**: Returns a vector of predictions `t̂ₛ` (`t̂ₛ` := "predictions of `t̄` based
    on source embedding `S̄`"), where `t̂ₛ[i]` is the prediction for `t̄[i]`. It assumes
    pre-embedded data which have been correctly time-aligned using a joint embedding
    (see [`embed`](@ref)), i.e. such that `t̄[i]` and `S̄[i]` correspond to the same time
    index.
- **Second method**: Jointly embeds the `target` and `source` time series (according to
    `measure`) to obtain time-index aligned target timeseries `t̄` and source embedding
    `S̄` (which is now a [`Dataset`](@ref)).
    Then calls `predict(measure, t̄, S̄)` (the first method), and returns both the
    predictions `t̂ₛ`, observations `t̄` and their correspondence `ρ` according to `measure`.

## Description

For each `i ∈ {1, 2, …, N}` where `N = length(t) == length(s)`, we make the prediction
`t̂[i]` (an estimate of `t[i]`) based on a linear combination of `D + 1` other points in
`t`, where the selection of points and weights for the linear combination are determined
by the `D+1` nearest neighbors of the point `S̄[i]`. The details of point selection and
weights depend on `measure`.

*Note: Some [`CrossmapMeasure`](@ref)s may define more general mapping
procedures. If so, the algorithm is described in their docstring*.
"""
function predict(measure::CrossmapMeasure, t::AbstractVector, S̄::AbstractDataset)
    @assert length(S̄) == length(t)
    (; d, τ, w) = measure
    # Tree structure must be created for every L, because we can't include data
    # outside the considered time range.

    # The number of neighbors depend on the type of cross map measure. We could make
    # this a tunable parameter, but for now, just stick with dim(embedding) + 1.
    nnd = n_neighbors_simplex(measure) #dimension(S̄) + 1

    # TODO: maybe construct the tree at a higher level, and use a more elaborate skip
    # function? Not sure what is fastest. Need to experiment...
    tree = KDTree(S̄, Euclidean())
    nnidxs, ds = bulksearch(tree, S̄, NeighborNumber(nnd), Theiler(w + 1))

    t̂ₛ = zeros(length(S̄))
    u = zeros(MVector{nnd}) # one extra neighbor due to the extra coordinate from t
    w = zeros(MVector{nnd}) # one extra neighbor due to the extra coordinate from t
    for (i, (nnidxsᵢ, dᵢ)) in enumerate(zip(nnidxs, ds))
        # There are two cases where the distance(s) from the query point to its closest
        # neighbor(s) is indistinguishable from zero:
        #   1) when sampling with replacement, or
        #   2) when input data points are very similar.
        # When either is the case, we encounter division by zero when computing weights.
        # To circumvent zero-division, we simply add some artifical small distance
        # to all distances, so that there are no zero-distance neighbors.
        if !(first(dᵢ) > 0.0)
            for i = 1:nnd
                # One order of magnitude higher than smallest possible float
                dᵢ[i] += eps()*10
            end
        end
        u .= exp.(-dᵢ ./ dᵢ[1])
        w .= u ./ sum(u)
        # Predict using weights computed from source `s` applied to values from target `t`
        t̂ₛ[i] = sum(w .* t[nnidxsᵢ])
    end
    return t̂ₛ
end

"""
    CrossmapEnsemble(; measure::CrossmapMeasure, est::CrossmapEstimator, nreps::Int = 100)

A directive to compute an ensemble cross-map analysis, where `measure` (e.g.
[`ConvergentCrossMapping`](@ref)) is computed
using the given estimator `est` (e.g. [`RandomVectors`](@ref))
"""
Base.@kwdef struct CrossmapEnsemble{M, E} # todo: perhaps just use a more general Ensemble?
    measure::M
    est::E
    nreps::Int = 100
    function CrossmapEnsemble(measure::M, est::E; nreps = 100) where {M, E}
        new{M, E}(measure, est, nreps)
    end
end

include("compat.jl") # remove this file when 2.0 is released.
include("estimators/estimators.jl")
include("measures/ccm-like/measures.jl")
