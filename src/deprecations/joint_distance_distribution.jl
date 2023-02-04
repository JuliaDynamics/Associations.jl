
"""
    jdd(measure::JointDistanceDistribution), source, target)

Compute the joint distance distribution (Amigó & Hirata, 2018[^Amigo2018]) from `source`
to `target` using the provided `distance_metric`, with `B` controlling the number of
subintervals, `D` the embedding dimension and `τ` the embedding lag.

Amigó & Hirata uses a one-sample t-test to determine whether the resulting
distance distribution is skewed towards positive values. You can obtain this
behaviour by doing

```julia
using CausalityTools
x, y = rand(100), rand(100)
measure = JointDistanceDistribution(; D = 15)
Δjdd = jdd(measure, x, y)
OneSampleTTest(Δjdd, measure.μ)
````

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
        "Call `OneSampleTTest manually on the result instead with the desired `μ`, i.e.\n" *
        "  measure = JointDistanceDistribution(; μ = 0.0)\n" *
        "  res = jdd(measure, x, y)\n"*
        "  OneSampleTTest(res, measure.μ, tail = :right) # at 95% confidence level"
    )
    measure = JointDistanceDistribution(; kwargs...)
    Δjdd = jdd(measure, x, y)
    return OneSampleTTest(Δjdd, measure.μ)
end
