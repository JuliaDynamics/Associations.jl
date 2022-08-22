using Statistics
include("ccc_estimators.jl")

function icc_i(Xⁱ::AbstractVector{J}, X⁻::AbstractVector{J}, Y⁻::AbstractVector{J},
    estimator::CompressionComplexityCausalityEstimator) where {J <: Integer}
    return dc(Xⁱ, X⁻, estimator.algorithm) - dc(Xⁱ, X⁻, Y⁻, estimator.algorithm)
end

function icc_on_segments(x::AbstractVector{J}, y::AbstractVector{J},
        estimator::CompressionComplexityCausalityEstimator) where {J <: Integer}
    w, L, step = estimator.w, estimator.L, estimator.step
    windows = get_windows(x, w + L, step)
    segment_results = zeros(eltype(zero(1.0)), length(windows))
    for (i, window) in enumerate(windows)
        current_indices = window[estimator.w + 1]:window[estimator.w + estimator.L]
        past_indices = window[1:estimator.w]
        Xⁱ = @views x[current_indices]
        X⁻ = @views x[past_indices]
        Y⁻ = @views y[past_indices]
        segment_results[i] = dynamical_complexity(Xⁱ, X⁻, estimator.algorithm) -
            dynamical_complexity(Xⁱ, X⁻, Y⁻, estimator.algorithm)
    end
    return segment_results
end

function icc(x, y, estimator)
    return Statistics.mean(icc_on_segments(x, y, estimator))
end
