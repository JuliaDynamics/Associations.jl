include("ccc_estimators.jl")

function icc_i(Xⁱ::AbstractVector{J}, X⁻::AbstractVector{J}, Y⁻::AbstractVector{J}, 
    estimator::CompressionComplexityCausalityEstimator) where {J <: Integer}
    return dc(Xⁱ, X⁻, estimator.algorithm) - dc(Xⁱ, X⁻, Y⁻, estimator.algorithm)
end

function icc_on_segments(x::AbstractVector{J}, y::AbstractVector{J}, estimator::CompressionComplexityCausalityEstimator) where {J <: Integer}
    window_size = estimator.w + estimator.L
    windows = get_windows(x, window_size, estimator.step)
    segment_results = zeros(Float64, length(windows))
    for (i, window) in enumerate(windows)
        current_indices = window[est.w+1]:window[est.w+est.L]
        past_indices = window[1:est.w]
        Xⁱ = @views x[current_indices]
        X⁻ = @views x[past_indices]
        Y⁻ = @views y[past_indices]
        segment_results[i] = dc(Xⁱ, X⁻, estimator.algorithm) - 
            dc(Xⁱ, X⁻, Y⁻, estimator.algorithm)
    end
    return segment_results
end

function icc(x::AbstractVector{J}, y::AbstractVector{J}, 
        estimator::CompressionComplexityCausalityEstimator) where {J <: Integer}
    return mean(icc_on_segments(x, y, estimator))
end