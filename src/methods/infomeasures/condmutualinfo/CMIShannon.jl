using Accessors

export CMIShannon
import ComplexityMeasures: log_with_base

"""
    CMIShannon <: ConditionalMutualInformation
    CMIShannon(; base = 2)

The Shannon conditional mutual information (CMI) ``I^S(X; Y | Z)``.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence.
- Use with [`condmutualinfo`](@ref) to compute the raw conditional mutual information.

## Supported definitions

Consider random variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}``, given ``Z \\in \\mathbb{R}^{d_Z}``. The Shannon
conditional mutual information is defined as

```math
\\begin{align*}
I(X; Y | Z)
&= H^S(X, Z) + H^S(Y, z) - H^S(X, Y, Z) - H^S(Z) \\\\
&= I^S(X; Y, Z) + I^S(X; Y)
\\end{align*},
```

where ``I^S(\\cdot; \\cdot)`` is the Shannon mutual information [`MIShannon`](@ref),
and ``H^S(\\cdot)`` is the [`Shannon`](@ref) entropy.

Differential Shannon CMI is obtained by replacing the entropies by
differential entropies.

See also: [`condmutualinfo`](@ref).
"""
struct CMIShannon{E <: Shannon} <: ConditionalMutualInformation{E}
    e::E
    function CMIShannon(; base::T = 2) where {T <: Real}
        e = Shannon(; base)
        new{typeof(e)}(e)
    end
    function CMIShannon(e::E) where E <: Shannon
        new{E}(e)
    end
end

function estimate(measure::CMIShannon, est::OutcomeSpace, x, y, z)
    HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h(measure, est, x, y, z)
    return HXZ + HYZ - HXYZ - HZ
end


function estimate(measure::CMIShannon, est::DifferentialInfoEstimator, x, y, z)
    # Due to inconsistent API in ComplexityMeasures.jl, we have to treat
    # DifferentialInfoEstimator here. Because all measures in this package
    # have their own `base` field, it will conflict with `est.base` for
    # `DifferentialInfoEstimator`s. In these cases, we use `measure.base`,
    # and override the estimator base, by simply creating a copy of the
    # estimator with one field modified.
    if est isa DifferentialInfoEstimator && :base in fieldnames(typeof(est))
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

function estimate(measure::CMIShannon, est::Contingency{<:OutcomeSpace}, x...)
    return estimate(measure, contingency_matrix(est.est, x...))
end

function estimate(measure::CMIShannon, est::Contingency{<:Nothing}, x...)
    return estimate(measure, contingency_matrix(x...))
end

function estimate(
        measure::CMIShannon,
        pxyz::ContingencyMatrix{T, 3}) where T
    e = measure.e
    dx, dy, dz = size(pxyz)
    pxz = probabilities(pxyz, dims = [1, 3])
    pyz = probabilities(pxyz, dims = [2, 3])
    pz = probabilities(pxyz, dims = 3)
    cmi = 0.0
    log0 = log_with_base(e.base)
    for k in 1:dz
        pzₖ = pz[k]
        for j in 1:dy
            pyⱼzₖ = pyz[j, k]
            pyⱼzₖ > 0 || continue # leads to NaN
            for i in 1:dx
                pxᵢzₖ = pxz[i, k]
                pxᵢzₖ > 0 || continue # leads to NaN
                pxᵢyⱼzₖ = pxyz[i, j, k]
                inner = (pzₖ * pxᵢyⱼzₖ) / (pxᵢzₖ * pyⱼzₖ)
                if inner != 0.0
                    cmi += pxᵢyⱼzₖ * log0(inner)
                end
            end
        end
    end
    return cmi
end
