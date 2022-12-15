export MITsallis

"""
    MITsallis <: MutualInformation
    MITsallis(; base = 2, q = 1.5)

The [`Tsallis`](@ref) mutual information measure.

## Supported definitions

- [`TsallisH3`](@ref).

If used with a [`MutationalInformationEstimator`](@ref), then the estimator dictates the
definition (i.e. which formula is computed).
"""
struct MITsallis{E <: Tsallis} <: MutualInformation
    e::E
    function MITsallis(; q = 1.5, base = 2)
        e = Tsallis(; q, base)
        new{typeof(e)}(e)
    end
end

# Define default behaviour.
estimate(measure::MITsallis, est, x, y) = estimate(TsallisH3(), measure, est, x, y)
