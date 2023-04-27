using Accessors

export CMIRenyiJizba
"""
    CMIRenyiJizba <: ConditionalMutualInformation

The Rényi conditional mutual information ``I_q^{R_{J}}(X; Y | Z`` defined in Jizba et
al. (2012)[^Jizba2012].

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence.
- Use with [`condmutualinfo`](@ref) to compute the raw conditional mutual information.

## Definition

```math
I_q^{R_{J}}(X; Y | Z) = I_q^{R_{J}}(X; Y, Z) - I_q^{R_{J}}(X; Z),
```

where ``I_q^{R_{J}}(X; Z)`` is the [`MIRenyiJizba`](@ref) mutual information.

[^Jizba2012]:
    Jizba, P., Kleinert, H., & Shefaat, M. (2012). Rényi’s information transfer between
    financial time series. Physica A: Statistical Mechanics and its Applications,
    391(10), 2971-2989.
"""
struct CMIRenyiJizba{E <: Renyi} <: ConditionalMutualInformation{E}
    e::E
    function CMIRenyiJizba(; base = 2, q = 1.5)
        e = Renyi(; base, q)
        new{typeof(e)}(e)
    end
    function CMIRenyiJizba(e::E) where E <: Renyi
        new{E}(e)
    end
end

function estimate(measure::CMIRenyiJizba, est::Contingency, x, y, z)
    c = _contingency_matrix(measure, est, x, y, z)
    pxz = probabilities(c, dims = [1, 3])
    pyz = probabilities(c, dims = [2, 3])
    pz = probabilities(c, dims = 3)
    pxyz = probabilities(c)
    e = measure.e
    return entropy(e, pxz) + entropy(e, pyz) - entropy(e, pz) - entropy(e, pxyz)
end

function _contingency_matrix(measure::CMIRenyiJizba,
        est::Contingency{<:ProbabilitiesEstimator}, x, y, z)
    return contingency_matrix(est.est, x, y, z)
end
function _contingency_matrix(measure::CMIRenyiJizba, est::Contingency{<:Nothing}, x, y, z)
    return contingency_matrix(x, y, z)
end

function estimate(measure::CMIRenyiJizba, est::ProbabilitiesEstimator, x, y, z)
    HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h(measure, est, x, y, z)
    return HXZ + HYZ - HXYZ - HZ
end

function estimate(measure::CMIRenyiJizba, est::DifferentialEntropyEstimator, x, y, z)
    # Due to inconsistent API in ComplexityMeasures.jl, we have to treat
    # DifferentialEntropyEstimator here. Because all measures in this package
    # have their own `base` field, it will conflict with `est.base` for
    # `DifferentialEntropyEstimator`s. In these cases, we use `measure.base`,
    # and override the estimator base, by simply creating a copy of the
    # estimator with one field modified.
    if est isa DifferentialEntropyEstimator && :base in fieldnames(typeof(est))
        if est.base != measure.e.base
            mb = measure.e.base
            eb = est.base
            modified_est = Accessors.@set est.base = measure.e.base
            HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h(measure, modified_est, x, y, z)
        else
            HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h(measure, est, x, y, z)
        end
    else
        HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h(measure, est, x, y, z)
    end
    cmi = HXZ + HYZ - HXYZ - HZ
    return cmi
end
