"""
    MIRenyiJizba <: ConditionalMutualInformation

The Rényi mutual information ``I_q^{R_{J}}(X; Y)`` defined in
Jizba et al. (2012)[^Jizba2012].

## Definition

```math
I_q^{R_{J}}(X; Y) = S_q^{R}(X) + S_q^{R}(Y) - S_q^{R}(X, Y),
```

where ``S_q^{R}(\\cdot)`` and ``S_q^{R}(\\cdot, \\cdot)`` the [`Rényi`](@ref) entropy and
the joint Rényi entropy.

[^Jizba2012]:
    Jizba, P., Kleinert, H., & Shefaat, M. (2012). Rényi’s information transfer between
    financial time series. Physica A: Statistical Mechanics and its Applications,
    391(10), 2971-2989.
"""
struct MIRenyiJizba{E <: Renyi} <: MutualInformation{E}
    e::E
    function MIRenyiJizba(; q = 1.5, base = 2)
        e = Renyi(; q, base)
        new{typeof(e)}(e)
    end
    function MIRenyiJizba(e::E) where E <: Renyi
        new{E}(e)
    end
end

function estimate(measure::MIRenyiJizba, est::ProbOrDiffEst, x, y)
    HX, HY, HXY = marginal_entropies_mi3h(measure, est, x, y)
    return HX + HY - HXY
end
