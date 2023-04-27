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
function crossmap(measure::CCMLike, est::CrossmapEstimator, target, source)
    return last.(predict(measure, est, target, source))
end
function crossmap(measure::CCMLike, est::CrossmapEstimator{<:Integer}, target, source)
    return last(predict(measure, est, target, source))
end

function predict(measure::CCMLike, est::CrossmapEstimator, target::AbstractVector, source::AbstractVector)
    emb, idx_t̄, idxs_S̄ = embed(measure, target, source)
    S̄ = emb[:, idxs_S̄]
    t̄ = emb[:, idx_t̄]
    return predict(measure, est, t̄, S̄)
end

# The following methods assume pre-embedded data.
# =========================================================================================
function predict(measure::CCMLike, est::CrossmapEstimator, target::AbstractVector, source::AbstractStateSpaceSet)
    # Ensure equal-length input
    input_check(measure, target, source)

    n_libraries = length(est.libsizes)
    ρs = Vector{Tuple{Vector{<:Real}, <:Real}}(undef, n_libraries)
    for i = 1:n_libraries
        # Randomly or deterministically determined indices for the library points.
        inds = library_indices(measure, est, i, target, source)
        # Predict on the library (i.e. on selected subset of points).
        ρs[i] = subset_predict(measure, target, source, inds)
    end
    return ρs
end

function predict(measure::CCMLike, est::CrossmapEstimator{<:Integer}, target::AbstractVector, source::AbstractStateSpaceSet)
    # Ensure equal-length input
    input_check(measure, target, source)
    inds = library_indices(measure, est, 1, target, source)
    ρ = subset_predict(measure, target, source, inds)
    return ρ
end

function subset_predict(measure::CCMLike, target, source, inds)
    S̄ = @views source[inds]
    t̄ = @views target[inds]
    t̂ₛ = predict(measure, t̄, S̄)
    ρ = measure.f(t̄, t̂ₛ)
    return t̂ₛ, ρ
end

"""
    library_indices(measure::CCMLike, est::CrossmapEstimator, i::Int,  target, source)

Produce (randomly, if relevant) the `i`-th subset of indices for a `CrossmapEstimator`
that is being applied to `target` and `source`.
"""
function library_indices end

function library_indices(measure::CCMLike, est::RandomVectors, i::Int, target, args...)
    N = length(target)
    L = est.libsizes[i]
    return library_indices(measure, est, N, L)
end
function library_indices(measure::CCMLike, est::RandomVectors, N::Int, L::Int)
    return sample(est.rng, 1:N, L; replace = est.replace)
end

function library_indices(measure::CCMLike, est::RandomSegment, i::Int, target, args...)
    N = length(target)
    L = est.libsizes[i]
    Lmax = max_segmentlength(measure, target)
    L <= Lmax ||
        throw(ArgumentError("L ($L) > Lmax ($Lmax). Use a smaller segment length (some points are lost when embedding)."))
    library_indices(measure, est, N, L)
end
function library_indices(measure::CCMLike, est::RandomSegment, N::Int, L::Int)
    startidx = sample(est.rng, 1:(N - L)) # random segment starting point
    return startidx:startidx+L-1
end

function library_indices(measure::CCMLike, est::ExpandingSegment, i::Int, target, args...)
    Lmax = max_segmentlength(measure, target)
    L = est.libsizes[i]
    L <= Lmax || throw(ArgumentError("L ($L) > Lmax ($Lmax). Use a smaller segment length (some points are lost when embedding)."))
    return library_indices(measure, est, length(target), L)
end
function library_indices(measure::CCMLike, est::ExpandingSegment, N::Int, L::Int)
    return 1:L
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
function crossmap(ensemble::Ensemble{<:CCMLike, <:CrossmapEstimator{<:Integer, R}},
        target::AbstractVector, source::AbstractVector) where R
    (; measure, est, nreps) = ensemble
    emb, idx_t̄, idxs_S̄ = embed(measure, target, source)
    S̄ = emb[:, idxs_S̄]
    t̄ = emb[:, idx_t̄]

    ρs = zeros(nreps)
    for i = 1:nreps
        inds = library_indices(measure, est, 1, t̄, S̄)
        ρ = last(subset_predict(measure, target, S̄, inds))
        ρs[i] = ρ
    end
    return ρs
end

function crossmap(ensemble::Ensemble{<:CCMLike, <:CrossmapEstimator},
        target::AbstractVector, source::AbstractVector)
    (; measure, est, nreps) = ensemble
    libsizes = est.libsizes
    emb, idx_t̄, idxs_S̄ = embed(measure, target, source)
    S̄ = emb[:, idxs_S̄]
    t̄ = emb[:, idx_t̄]
    N = length(t̄)

    ρs = [zeros(nreps) for j in eachindex(libsizes)]
    for (j, L) in enumerate(libsizes)
        for i = 1:nreps
            inds = library_indices(measure, est, N, L)
            ρs[j][i] = last(subset_predict(measure, target, S̄, inds))
        end
    end
    return ρs
end
