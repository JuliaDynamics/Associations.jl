
"""
    jdd(measure::JointDistanceDistribution, source, target) → Δ

Compute the joint distance distribution (Amigó & Hirata, 2018[^Amigo2018]) from `source`
to `target` using the given [`JointDistanceDistribution`](@ref) measure.

Returns the distribution `Δ` from the paper directly ([example](@ref quickstart_jdd)).
Use [`JointDistanceDistributionTest`](@ref) to perform a formal indepencence test.

[^Amigo2018]:
    Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from multivariate
    flows by the joint distance distribution." Chaos: An Interdisciplinary Journal of
    Nonlinear Science 28.7 (2018): 075302.
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
