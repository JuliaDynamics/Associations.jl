export CCCEstimator, CompressionComplexityCausalityEstimator

"""
    CompressionComplexityCausalityEstimator(
        algorithm::CompressionComplexityAlgorithm = EffortToCompress(),
        w::In = 15,
        L::Int = 30)

A `CompressionComplexityCausalityEstimator` (`CCCEstimator` for short) stores 
parameters used to calculate compute interventional complexity causality:
which compression `algorithm` to use, as well as the window sizes `w` 
(for present values) and `L` (for past values) of the relevant variables.

## Input data requirements

All estimators assume that the input time series are *pre-symbolized/binned* integer-valued 
time series. Thus, the parameter `B` (number of bins) from Kathpalia & Nararaj 
(2019)[^Kathpalia2019] is not present here.

See also: [`icc`](@ref), [`EffortToCompress`](@ref).

[^Kathpalia2019]: Kathpalia, A., & Nagaraj, N. (2019). Data-based intervention approach for Complexity-Causality measure. PeerJ Computer Science, 5, e196.
"""
Base.@kwdef struct CompressionComplexityCausalityEstimator{ALG} <: InterventionalComplexityCausalityEstimator
    """
    The algorithm used to compute the compression complexity, e.g. an instance of 
    [`EffortToCompress`](@ref).
    """
    algorithm::ALG = EffortToCompress(normalize = true)

    """ 
    The block size for `Xⁱ`, which is a block of current values of time series `X`. The 
    default value is set to 15, as recommended in Kathpalia & Nagaraj (2019).
    """
    w = 15

    """ 
    The block size for `X⁻`, `Y⁻`, `Z⁻`, ... , which are blocks of past values of 
    `X`, `Y`, `Z`, and so on, whose values are taken from the immediate past of the `Xⁱ`
    block. The default value is `L = 30`, but for real applications, this value 
    should determined as described in table S2 in the supplement of Kathpalia & Nagaraj 
    (2019).
    """
    L = 30

    step = 1
end
const CCCEstimator = CompressionComplexityCausalityEstimator

