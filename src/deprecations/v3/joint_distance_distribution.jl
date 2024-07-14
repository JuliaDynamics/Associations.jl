export jdd 

"""
    jdd(measure::JointDistanceDistribution, source, target) → Δ

Compute the joint distance distribution [Amigo2018](@cite) from
`source` to `target` using the given [`JointDistanceDistribution`](@ref) measure.

Returns the distribution `Δ` from the paper directly ([example](@ref quickstart_jdd)).
Use [`JointDistanceDistributionTest`](@ref) to perform a formal indepencence test.
"""
function jdd(source, target; kw...)
    @warn(
        "Convenience function `jdd` is deprecated. " *
        "Use `association(JointDistanceDistribution(; kwargs...), source, target) instead."
    )
    association(JointDistanceDistribution(; kw...), source, target)
end

function jdd(measure::JointDistanceDistribution, source, target)
    @warn(
        "Convenience function `jdd` is deprecated. " *
        "Use `association(JointDistanceDistribution(; kwargs...), source, target) instead."
    )
    return association(measure, source, target)
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
