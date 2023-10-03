using ComplexityMeasures: Tsallis

export MITsallisMartin

"""
    MITsallisMartin <: BivariateInformationMeasure
    MITsallisMartin(; base = 2, q = 1.5)

The discrete Tsallis mutual information from [Martin2004](@citet).

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise
    dependence.
- Use with [`mutualinfo`](@ref) to compute the raw mutual information.

## Description

Martin et al.'s Tsallis mutual information between variables ``X \\in \\mathbb{R}^{d_X}``
and ``Y \\in \\mathbb{R}^{d_Y}`` is defined as

```math
I_{\\text{Martin}}^T(X, Y, q) := H_q^T(X) + H_q^T(Y) - (1 - q) H_q^T(X) H_q^T(Y) - H_q(X, Y),
```

where ``H^S(\\cdot)`` and ``H^S(\\cdot, \\cdot)`` are the marginal and joint Shannon
entropies, and `q` is the [`Tsallis`](@ref)-parameter.

See also: [`mutualinfo`](@ref).
"""
struct MITsallisMartin{E <: Tsallis} <: MutualInformation
    e::E
    function MITsallisMartin(; q = 1.5, base = 2)
        e = Tsallis(; q, base)
        new{typeof(e)}(e)
    end
end

# This is definition 3 in Martin et al. (2004), but with pᵢ replaced by the joint
# distribution and qᵢ replaced by the product of the marginal distributions.
function information(definition::MITsallisMartin, pxy::Probabilities{T, 2}) where T
    e = definition.e
    q = definition.e.q
    q != 1 || throw(ArgumentError("`MITsallisMartin` for q=$(q) not defined."))
    px = marginal(pxy, dims = 1)
    py = marginal(pxy, dims = 2)

    mi = 0.0
    for (i, pxᵢ) in enumerate(px.p)
        for (j, pyⱼ) in enumerate(py.p)
            pxyᵢⱼ = pxy[i, j]
            mi += pxyᵢⱼ^q / (pxᵢ^(q - 1) * pyⱼ^(q - 1))
        end
    end
    f = 1 / (q - 1)
    return f * (1 - mi)
end