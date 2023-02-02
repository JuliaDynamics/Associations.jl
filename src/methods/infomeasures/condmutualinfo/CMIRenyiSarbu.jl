export CMIRenyiSarbu

"""
    CMIRenyiSarbu <: ConditionalMutualInformation
    CMIRenyiSarbu(; base = 2, definition = CMIRenyiSarbuSarbu())

The Rényi conditional mutual information from Sarbu (2014)[^Sarbu2014]).

## Discrete description

Assume we observe three discrete random variables ``X``, ``Y`` and ``Z``.
Sarbu (2014) defines discrete conditional Rényi mutual information as the conditional
Rényi ``\\alpha``-divergence between the conditional joint probability mass function
``p(x, y | z)`` and the product of the conditional marginals, ``p(x |z) \\cdot p(y|z)``:

```math
I(X, Y; Z)^R_q =
\\dfrac{1}{q-1} \\sum_{z \\in Z} p(Z = z)
\\log \\left(
    \\sum{x \\in X}\\sum{y \\in Y}
    \\dfrac{p(x, y|z)^q}{\\left( p(x|z)\\cdot p(y|z) \\right)^{q-1}}
\\right)
```

[^Sarbu2014]: Sarbu, S. (2014, May). Rényi information transfer: Partial Rényi transfer
    entropy and partial Rényi mutual information. In 2014 IEEE International Conference
    on Acoustics, Speech and Signal Processing (ICASSP) (pp. 5666-5670). IEEE.

See also: [`condmutualinfo`](@ref).
"""
struct CMIRenyiSarbu{E <: Renyi} <: ConditionalMutualInformation{E}
    e::E
    function CMIRenyiSarbu(; base = 2, q = 1.5)
        e = Renyi(; base, q)
        new{typeof(e)}(e)
    end
end

function estimate(measure::CMIRenyiSarbu, est::Contingency{<:ProbabilitiesEstimator}, x...)
    return estimate(measure, contingency_matrix(est.est, x...))
end

function estimate(measure::CMIRenyiSarbu, est::Contingency{<:Nothing}, x...)
    return estimate(measure, contingency_matrix(x...))
end

function estimate(
        measure::CMIRenyiSarbu,
        pxyz::ContingencyMatrix{T, 3}) where T
    e = measure.e
    q = e.q
    dx, dy, dz = size(pxyz)
    pxz = probabilities(pxyz, dims = [1, 3])
    pyz = probabilities(pxyz, dims = [2, 3])
    pz = probabilities(pxyz, dims = 3)
    cmi = 0.0
    logb = log_with_base(e.base)
    for k in 1:dz
        pzₖ = pz[k]
        inner = 0.0
        pzₖ != 0.0 || continue
        for j in 1:dy
            pyⱼzₖ = pyz[j, k]
            for i in 1:dx
                pxᵢzₖ = pxz[i, k]
                pxᵢyⱼzₖ = pxyz[i, j, k]
                f = (pxᵢzₖ / pzₖ) * (pyⱼzₖ / pzₖ)
                if f != 0.0 && !isnan(f)
                    inner += (pxᵢyⱼzₖ / pzₖ)^q / f^(q-1)
                end
            end
        end
        cmi += pzₖ * logb(inner)
    end
    return 1 / (1 - q) * cmi
end
