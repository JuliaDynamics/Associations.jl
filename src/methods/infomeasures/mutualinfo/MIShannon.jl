export MIShannon

"""
    MIShannon <: MutualInformation
    MIShannon(; base = 2)

The Shannon mutual information ``I^S(X; Y)``.

## Direct definition: discrete

Assume we observe samples
``\\bar{\\bf{X}}_{1:N_y} = \\{\\bar{\\bf{X}}_1, \\ldots, \\bar{\\bf{X}}_n \\}`` and
``\\bar{\\bf{Y}}_{1:N_x} = \\{\\bar{\\bf{Y}}_1, \\ldots, \\bar{\\bf{Y}}_n \\}`` from
two discrete random variables ``X`` and ``Y`` with finite supports
``\\mathcal{X} = \\{ x_1, x_2, \\ldots, x_{M_x} \\}`` and
``\\mathcal{Y} = y_1, y_2, \\ldots, x_{M_y}``.

In practice, we estimate the double sum
``I(X; Y) = \\sum_{x \\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y) log \\left( \\dfrac{p(x, y)}{p(x)p(y)} \\right)``
by replacing the true probabilities ``p(x_i)``, ``p(y_j)``, and ``p(x_i, y_j)`` by their
maximum likelihood empirical estimates ``\\hat{p}(x_i) = \\frac{n(x_i)}{N_x}``,
``\\hat{p}(y_i) = \\frac{n(y_j)}{N_y}``,
and ``\\hat{p}(x_i, x_j) = \\frac{n(x_i)}{N}``, where ``N = N_x N_y``.


```math
\\begin{align*}
\\hat{I}_{DS}(X; Y) &=
 \\sum_{x_i \\in \\mathcal{X}, y_i \\in \\mathcal{Y}} p(x_i, y_j) \\log \\left( \\dfrac{p(x_i, y_i)}{p(x_i)p(y_j)} \\right) \\\\
&= p(x_1, y_1) \\log p(x_1, y_1) + p(x_1, y_2) \\log p(x_1, y_2) + \\cdots + p(x_1, y_{M_y}) \\log p(x_1, y_{M_y}) \\\\
&+ p(x_2, y_1) \\log p(x_2, y_1) + p(x_2, y_2) \\log p(x_2, y_2) + \\cdots + p(x_2, y_{M_y}) \\log p(x_2, y_{M_y}) \\\\
&+ \\cdots \\\\
&+ p(x_{M_x}, y_1) \\log p(x_{M_x}, y_1) + p(x_{M_x}, y_2) \\log p(x_{M_x}, y_2) + \\cdots + p(x_{M_x}, y_1) \\log p(x_{M_x}, y_{M_y}) \\\\
\\end{align*}
```

This definition is used by [`mutualinfo`](@ref) when called with a
[`ContingencyMatrix`](@ref).

## H3 definition: discrete and continuous

The H3-definition and equivalent definition of Shannon mutual information is

- Continuous case: ``I^S(X; Y) = H^S(X) + H_q^S(Y) - H^S(X, Y)``
- Discrete case: ``I^S(X; Y) = H^S(X) + H_q^S(Y) - H^S(X, Y)``,

Here, ``H^S(\\cdot)`` and ``H^S(\\cdot, \\cdot)`` are the marginal and joint discrete
Shannon entropies, and ``h^S(\\cdot)`` and ``h^S(\\cdot, \\cdot)`` are the corresponding
differential entropies. To compute ``I^S(X; Y)``, we simply compute each marginal
entropy separately and sum them according to the formuals above. The discrete
variant is equivalent to the direct formulation above.

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

function estimate(measure::MIShannon, pxy::ContingencyMatrix{T, 2}) where {T}
    e = measure.e
    px = probabilities(pxy, 1)
    py = probabilities(pxy, 2)
    mi = 0.0
    log0 = log_with_base(e.base)
    for i in eachindex(px.p)
        pxᵢ = px[i]
        for j in eachindex(py.p)
            pyⱼ = py[j]
            pxyᵢⱼ = pxy[i, j]
            mi += pxyᵢⱼ * log0(pxyᵢⱼ / (pxᵢ * pyⱼ))
        end
    end
    return mi
end

function estimate(measure::MIShannon, est::ProbOrDiffEst, x, y)
    HX, HY, HXY = marginal_entropies_mi3h(measure, est, x, y)
    return HX + HY - HXY
end
