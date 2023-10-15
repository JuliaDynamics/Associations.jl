using StatsBase: sample, sample!
using StateSpaceSets: dimension
using StateSpaceSets: AbstractStateSpaceSet

CCMLike = Union{ConvergentCrossMapping, PairwiseAsymmetricInference}

# -----------------------------------------------------------------------------------------
# Generic methods that operates on the entire input. Useful for reconstructing figures
# from e.g. Sugihara et al.
# -----------------------------------------------------------------------------------------
crossmap(measure::CCMLike, target, source) = last(predict(measure, target, source))
function predict(measure::CrossmapMeasure, target::AbstractVector, source::AbstractVector)
    emb, idx_t̄, idxs_S̄ = embed(measure, target, source)
    S̄ = emb[:, idxs_S̄]
    t̄ = emb[:, idx_t̄]
    t̄ₛ = predict(measure, t̄, S̄);
    return t̄ₛ, t̄, measure.f(t̄, t̄ₛ)
end

# -----------------------------------------------------------------------------------------
# Estimator-specific implementations. Assumed pre-embedded data.
# -----------------------------------------------------------------------------------------

# Wrappers for timeseries inputs that ensure embeddings are done correctly.
# =========================================================================================
function crossmap(est::CrossmapEstimator{<:CCMLike}, target, source)
    return last.(predict(est, target, source))
end
function crossmap(est::CrossmapEstimator{<:CCMLike, <:Integer}, target, source)
    return last(predict(est, target, source))
end

function predict(est::CrossmapEstimator{<:CCMLike}, target::AbstractVector, source::AbstractVector)
    emb, idx_t̄, idxs_S̄ = embed(est.definition, target, source)
    S̄ = emb[:, idxs_S̄]
    t̄ = emb[:, idx_t̄]
    return predict(est, t̄, S̄)
end

# The following methods assume pre-embedded data.
# =========================================================================================
function predict(est::CrossmapEstimator{<:CCMLike}, target::AbstractVector, source::AbstractStateSpaceSet)
    # Ensure equal-length input
    input_check(est.definition, target, source)

    n_libraries = length(est.libsizes)
    ρs = Vector{Tuple{Vector{<:Real}, <:Real}}(undef, n_libraries)
    for i = 1:n_libraries
        # Randomly or deterministically determined indices for the library points.
        inds = library_indices(est, i, target, source)
        # Predict on the library (i.e. on selected subset of points).
        ρs[i] = subset_predict(est, target, source, inds)
    end
    return ρs
end

function predict(
        est::CrossmapEstimator{<:CCMLike, <:Integer}, 
        target::AbstractVector, 
        source::AbstractStateSpaceSet
    )
    definition = est.definition

    # Ensure equal-length input
    input_check(definition, target, source)
    inds = library_indices(est, 1, target, source)
    ρ = subset_predict(est, target, source, inds)
    return ρ
end

function subset_predict(est::CrossmapEstimator{<:CCMLike}, target, source, inds)
    definition = est.definition

    S̄ = @views source[inds]
    t̄ = @views target[inds]
    t̂ₛ = predict(definition, t̄, S̄)
    ρ = definition.f(t̄, t̂ₛ)
    return t̂ₛ, ρ
end

"""
    library_indices(measure::CCMLike, est::CrossmapEstimator, i::Int,  target, source)

Produce (randomly, if relevant) the `i`-th subset of indices for a `CrossmapEstimator`
that is being applied to `target` and `source`.
"""
function library_indices end

function library_indices(est::RandomVectors{<:CCMLike}, i::Int, target, args...)
    N = length(target)
    L = est.libsizes[i]
    return library_indices(est, N, L)
end
function library_indices(est::RandomVectors{<:CCMLike}, N::Int, L::Int)
    return sample(est.rng, 1:N, L; replace = est.replace)
end

function library_indices(est::RandomSegment{<:CCMLike}, i::Int, target, args...)
    N = length(target)
    L = est.libsizes[i]
    Lmax = max_segmentlength(est.definition, target)
    L <= Lmax ||
        throw(ArgumentError("L ($L) > Lmax ($Lmax). Use a smaller segment length (some points are lost when embedding)."))
    return library_indices(est, N, L)
end
function library_indices(est::RandomSegment{<:CCMLike}, N::Int, L::Int)
    startidx = sample(est.rng, 1:(N - L)) # random segment starting point
    return startidx:startidx+L-1
end

function input_check(measure::CCMLike, args...)
    ns = length.(args)
    all(ns .== maximum(ns)) || throw(ArgumentError("""\
        All inputs must have same lengths. \
        Use `embed` to ensure target time series and embedding are aligned.\
        """))
end

# # -----------------------------------------------------------------------------------------
# # Ensemble analysis. Repeats an analysis ensemble.nreps times. Takes care of the embedding.
# # -----------------------------------------------------------------------------------------
function crossmap(ensemble::Ensemble{<:CrossmapEstimator{<:CCMLike, <:Integer, R}},
        target::AbstractVector, source::AbstractVector) where R
    (; est, nreps) = ensemble
    definition = est.definition

    emb, idx_t̄, idxs_S̄ = embed(definition, target, source)
    S̄ = emb[:, idxs_S̄]
    t̄ = emb[:, idx_t̄]

    ρs = zeros(nreps)
    for i = 1:nreps
        inds = library_indices(est, 1, t̄, S̄)
        ρ = last(subset_predict(est, target, S̄, inds))
        ρs[i] = ρ
    end
    return ρs
end

function crossmap(
        ensemble::Ensemble{<:CrossmapEstimator},
        target::AbstractVector, 
        source::AbstractVector
    )
    (; est, nreps) = ensemble
    definition = est.definition

    libsizes = est.libsizes
    emb, idx_t̄, idxs_S̄ = embed(definition, target, source)
    S̄ = emb[:, idxs_S̄]
    t̄ = emb[:, idx_t̄]
    N = length(t̄)

    ρs = [zeros(nreps) for j in eachindex(libsizes)]
    for (j, L) in enumerate(libsizes)
        for i = 1:nreps
            inds = library_indices(est, N, L)
            ρs[j][i] = last(subset_predict(est, target, S̄, inds))
        end
    end
    return ρs
end
