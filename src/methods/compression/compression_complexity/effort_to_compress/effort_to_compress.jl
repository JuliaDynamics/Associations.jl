export EffortToCompress

"""
    EffortToCompress(; normalize = false)

The effort-to-compress (ETC; Nagaraj et al., 2013)[^Nagaraj2013] algorithm quantifies 
the compression complexity of a time series using non-sequential recursive pair 
substitution (NSRPS).

Used with [`complexity`](@ref) to compute ETC and with [`complexity_normalized`](@ref)
to compute normalized ETC. 
Input time series must be integer-valued. Real-valued or categorical time series 
must be symbolized (i.e. represented by integers) before computing the ETC. The number 
of symbols should be relatively low, i.e. `n_symbols << length(x)`, where `x` is the
input time series. If applied to two time series (e.g. 
`complexity(EffortToCompress(), x1, x2)`), then (bivariate) joint ETC (Kathpalia 
and Nagaraj, 2019)[^Kathpalia2019] is computed.

## Normalization

If `normalize == false`, then the quantity computed is the the number of compression steps 
it takes for the symbol sequence to reach zero entropy (either a constant sequence, 
or a length-1 sequence).

For a length-`N` sequence, the maximum number of possible compression steps is `N - 1`.
If `normalize == true`, then the ETC value is divided by `N - 1`, yielding a number on
`[0, 1]`, where `0` indicates minimal compression complexity and `1` indicates maximal
compression complexity.

## Description

The ETC algorithm parses a symbol sequence from left to right and finds the symbol with the
highest frequency of occurrence. *Note: In the case of multiple symbols attaining the 
maximum frequency, our implementation does not guarantee that the replaced symbol is 
the symbol which first occurs in the sequence. This behaviour differs from the examples in 
Nagaraj (2013) and Kathpalia (2019), where they seem to always replace the first-occurring
symbol. This does not affect the number of compression steps.*

[^Nagaraj2013]: Nagaraj, N., Balasubramanian, K., & Dey, S. (2013). A new complexity
    measure for time series analysis and classification. The European Physical Journal
    Special Topics, 222(3), 847-860.
[^Kathpalia2019]: Kathpalia, A., & Nagaraj, N. (2019). Data-based intervention approach
    for Complexity-Causality measure. PeerJ Computer Science, 5, e196.
"""
struct EffortToCompress <: CompressionComplexityEstimator end

include("etc_utils.jl")
include("etc_univariate.jl")
include("etc_multivariate.jl")
include("etc_joint.jl")
