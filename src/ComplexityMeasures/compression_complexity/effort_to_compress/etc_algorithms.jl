export EffortToCompress

"""
    EffortToCompress(; normalize = false)

The effort-to-compress (ETC; Nagaraj et al., 2013)[^Nagaraj2013] algorithm quantifies the compression complexity
of a time series.

If applied to two time series, the (bivariate) joint ETC as described in [^Kathpalia2019]
can also be computed.

## Normalization

If `normalize == false`, then computes the number of compression steps `N` it takes for the
symbol sequence to reach zero entropy (either a constant sequence, or a length-1 sequence).

For a length-`N` sequence, the maximum number of possible compression steps is `N-1`.
If `normalize == true`, then the ETC value is normalized to `N-1`, yielding a number on
`[0, 1]`, where `0` indicates minimal compression complexity and `1` indicates maximal
compression complexity.

[^Nagaraj2013]: Nagaraj, N., Balasubramanian, K., & Dey, S. (2013). A new complexity
    measure for time series analysis and classification. The European Physical Journal
    Special Topics, 222(3), 847-860.
[^Kathpalia2019]: Kathpalia, A., & Nagaraj, N. (2019). Data-based intervention approach
    for Complexity-Causality measure. PeerJ Computer Science, 5, e196.
"""
struct EffortToCompress <: CompressionComplexityAlgorithm
    normalize::Bool

    function EffortToCompress(; normalize = false)
        return new(normalize)
    end
end
