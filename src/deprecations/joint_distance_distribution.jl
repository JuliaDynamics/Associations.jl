
"""
    jdd(measure::JointDistanceDistribution, source, target) → Δ

Compute the joint distance distribution [Amigo2018](@cite) from
`source` to `target` using the given [`JointDistanceDistribution`](@ref) measure.

Returns the distribution `Δ` from the paper directly ([example](@ref quickstart_jdd)).
Use [`JointDistanceDistributionTest`](@ref) to perform a formal independence test.
"""
function jdd(source, target; kw...)
    if !isempty(kw)
        @warn(
            "Providing keywords to `jdd` is deprecated. " *
            "Use `jdd(JointDistanceDistribution(; kwargs...), source, target) instead of " *
            "`jdd(source, target; kwargs...)`"
        )
    end
    estimate(JointDistanceDistribution(; kw...), source, target)
end

function jdd(measure::JointDistanceDistribution, source, target)
    return estimate(measure, source, target)
end

function jdd(::Type{OneSampleTTest}, x, y; kwargs...)
    @warn(
        "jdd(::OneSampleTTest, x, y; kwargs...) is deprecated. " *
        "Instead, do\n" *
        "  measure = JointDistanceDistribution()\n" *
        "  independence(JointDistanceDistributionTest(measure), x, y)\n"
    )
    measure = JointDistanceDistribution(; kwargs...)
    Δjdd = jdd(measure, x, y)
    return OneSampleTTest(Δjdd, measure.μ)
end
