using ComplexityMeasures: log_with_base
using Accessors

export MIShannon

"""
    MIShannon <: MutualInformation
    MIShannon(; base = 2)

The Shannon mutual information ``I^S(X; Y)``.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence.
- Use with [`mutualinfo`](@ref) to compute the raw mutual information.

## Discrete definition

There are many equivalent formulations of discrete Shannon mutual information. In this
package, we currently use the double-sum and the three-entropies formulations.

### Double sum formulation

Assume we observe samples
``\\bar{\\bf{X}}_{1:N_y} = \\{\\bar{\\bf{X}}_1, \\ldots, \\bar{\\bf{X}}_n \\}`` and
``\\bar{\\bf{Y}}_{1:N_x} = \\{\\bar{\\bf{Y}}_1, \\ldots, \\bar{\\bf{Y}}_n \\}`` from
two discrete random variables ``X`` and ``Y`` with finite supports
``\\mathcal{X} = \\{ x_1, x_2, \\ldots, x_{M_x} \\}`` and
``\\mathcal{Y} = y_1, y_2, \\ldots, x_{M_y}``.
The double-sum estimate is obtained by replacing the double sum

```math
\\hat{I}_{DS}(X; Y) =
 \\sum_{x_i \\in \\mathcal{X}, y_i \\in \\mathcal{Y}} p(x_i, y_j) \\log \\left( \\dfrac{p(x_i, y_i)}{p(x_i)p(y_j)} \\right)
```

where  ``\\hat{p}(x_i) = \\frac{n(x_i)}{N_x}``, ``\\hat{p}(y_i) = \\frac{n(y_j)}{N_y}``,
and ``\\hat{p}(x_i, x_j) = \\frac{n(x_i)}{N}``, and ``N = N_x N_y``.
This definition is used by [`mutualinfo`](@ref) when called with a
[`ContingencyMatrix`](@ref).

### Three-entropies formulation

An equivalent formulation of discrete Shannon mutual information is

```math
I^S(X; Y) = H^S(X) + H_q^S(Y) - H^S(X, Y),
```

where ``H^S(\\cdot)`` and ``H^S(\\cdot, \\cdot)`` are the marginal and joint discrete
Shannon entropies. This definition is used by [`mutualinfo`](@ref) when called with a
[`ProbabilitiesEstimator`](@ref).

## Differential mutual information

One possible formulation of differential Shannon mutual information is

```math
I^S(X; Y) = h^S(X) + h_q^S(Y) - h^S(X, Y),
```

where ``h^S(\\cdot)`` and ``h^S(\\cdot, \\cdot)`` are the marginal and joint
differential Shannon entropies. This definition is used by [`mutualinfo`](@ref) when
called with a [`DifferentialInfoEstimator`](@ref).

See also: [`mutualinfo`](@ref).
"""
struct MIShannon{E <: Shannon} <: MutualInformation{E}
    e::E
    function MIShannon(; base::Real = 2)
        e = Shannon(; base)
        new{typeof(e)}(e)
    end
    function MIShannon(e::Shannon)
        new{typeof(e)}(e)
    end
end

function estimate(measure::MIShannon, est::Contingency{<:ProbabilitiesEstimator}, x...)
    return estimate(measure, contingency_matrix(est.est, x...))
end

function estimate(measure::MIShannon, est::Contingency{<:Nothing}, x...)
    return estimate(measure, contingency_matrix(x...))
end

function estimate(measure::MIShannon, pxy::ContingencyMatrix{T, 2}) where {T}
    e = measure.e
    px = probabilities(pxy, dims = 1)
    py = probabilities(pxy, dims = 2)
    mi = 0.0
    logb = log_with_base(e.base)
    for i in eachindex(px.p)
        pxᵢ = px[i]
        for j in eachindex(py.p)
            pyⱼ = py[j]
            pxyᵢⱼ = pxy[i, j]
            inner = pxyᵢⱼ / (pxᵢ * pyⱼ)
            if inner != 0.0
                mi += pxyᵢⱼ * logb(inner)
            end
        end
    end
    return mi
end

function estimate(measure::MIShannon, est::DiscreteOrDifferentialInfoEstimator, x, y)
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
            HX, HY, HXY = marginal_entropies_mi3h(measure, modified_est, x, y)
        else
            HX, HY, HXY = marginal_entropies_mi3h(measure, est, x, y)
        end
    else
        HX, HY, HXY = marginal_entropies_mi3h(measure, est, x, y)
    end
    mi = HX + HY - HXY
    return mi
end
