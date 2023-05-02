using StateSpaceSets
import ComplexityMeasures: complexity

function complexity(algorithm::EffortToCompress, x::AbstractStateSpaceSet{D, T}, 
        alphabet_size::Int) where {D, T}
    alphabet_size >= 2 || throw(ArgumentError("Alphabet size must be at least 2."))

    # Assuming that the alphabet size was the same when symbolizing all variables `xᵢ ∈ x`,
    # encode each unique state vector as a unique integer, so that the multivariate time
    # series becomes a univariate symbol time series.
    encoded_x = symbol_sequence(x, alphabet_size)
    return compression_complexity(encoded_x, algorithm)
end

# function compression_complexity(sw::ConstantWidthSlidingWindow{<:CompressionComplexityEstimator}, 
#         x::AbstractStateSpaceSet{D, T},
#         alphabet_size::Int) where {D, T}
#     alphabet_size >= 2 || throw(ArgumentError("Alphabet size must be at least 2."))

#     encoded_x = symbol_sequence(x, alphabet_size)
#     windows = get_windows(x, sw)
#     return @views [compression_complexity(encoded_x[w], sw.estimator) for w in windows]
# end
