using ComplexityMeasures: Renyi

export MIRenyiSarbu 

"""
    MIRenyiSarbu <: BivariateInformationMeasure
    MIRenyiSarbu(; base = 2, q = 1.5)

The discrete Rényi mutual information from [Sarbu2014](@citet).

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence.
- Use with [`mutualinfo`](@ref) to compute the raw mutual information.

## Description

Sarbu (2014) defines discrete Rényi mutual information as the
Rényi ``\\alpha``-divergence between the conditional joint probability mass function
``p(x, y)`` and the product of the conditional marginals, ``p(x) \\cdot p(y)``:

```math
I(X, Y)^R_q =
\\dfrac{1}{q-1}
\\log \\left(
    \\sum_{x \\in X, y \\in Y}
    \\dfrac{p(x, y)^q}{\\left( p(x)\\cdot p(y) \\right)^{q-1}}
\\right)
```
See also: [`mutualinfo`](@ref).
"""
struct MIRenyiSarbu{E <: Renyi} <: MutualInformation
    e::E
    function MIRenyiSarbu(; q = 1.5, base = 2)
        e = Renyi(; q, base)
        new{typeof(e)}(e)
    end
end

function information(definition::MIRenyiSarbu, pxy::Probabilities{T, 2}) where T
    px = marginal(pxy, dims = 1)
    py = marginal(pxy, dims = 2)
    e = definition.e
    q = e.q

    mi = 0.0
    for i in eachindex(px.p)
        for j in eachindex(py.p)
            pxyᵢⱼ = pxy[i, j]
            mi += pxyᵢⱼ^q / ((px[i] * py[j])^(q - 1))
        end
    end
    if mi == 0
        return 0.0
    else
        return _convert_logunit(1 / (q - 1) * log(mi), ℯ, e.base)
    end
end