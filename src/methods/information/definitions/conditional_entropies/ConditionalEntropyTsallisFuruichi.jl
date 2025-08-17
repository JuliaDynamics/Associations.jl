export ConditionalEntropyTsallisFuruichi

"""
    ConditionalEntropyTsallisFuruichi <: ConditionalEntropy
    ConditionalEntropyTsallisFuruichi(; base = 2, q = 1.5)

Furuichi (2006)'s discrete Tsallis conditional entropy definition.

## Usage 

- Use with [`association`](@ref) to compute the Tsallis-Furuichi conditional entropy between 
    two variables.

## Compatible estimators

- [`JointProbabilities`](@ref)

## Definition

Furuichi's Tsallis conditional entropy between discrete random variables
``X`` and ``Y`` with finite ranges ``\\mathcal{X}`` and ``\\mathcal{Y}`` is defined as

```math
H_q^T(X | Y) = -\\sum_{x \\in \\mathcal{X}, y \\in \\mathcal{Y}}
p(x, y)^q \\log_q(p(x | y)),
```

``\\ln_q(x) = \\frac{x^{1-q} - 1}{1 - q}`` and ``q \\neq 1``. For ``q = 1``, ``H_q^T(X | Y)`` reduces to the Shannon conditional
entropy:

```math
H_{q=1}^T(X | Y) = -\\sum_{x \\in \\mathcal{X}, y \\in \\mathcal{Y}} =
p(x, y) \\log(p(x | y))
```

If any of the entries of the marginal distribution for `Y` are zero, or the q-logarithm 
is undefined for a particular value, then the measure is undefined and `NaN` is returned.

## Estimation

- [Example 1](@ref example_ConditionalEntropyTsallisFuruichi_JointProbabilities_CodifyVariables_UniqueElements): 
    [`JointProbabilities`](@ref) estimator with[`CodifyVariables`](@ref) discretization and 
    [`UniqueElements`](@extref ComplexityMeasures) outcome space on categorical data.
- [Example 2](@ref example_ConditionalEntropyTsallisFuruichi_JointProbabilities_CodifyPoints_UniqueElementsEncoding): 
    [`JointProbabilities`](@ref) estimator with [`CodifyPoints`](@ref) discretization and [`UniqueElementsEncoding`](@ref)
    encoding of points on numerical data.
"""
Base.@kwdef struct ConditionalEntropyTsallisFuruichi{B,Q} <: ConditionalEntropy
    base::B = 2
    q::Q = 1.5
end

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:ConditionalEntropyTsallisFuruichi}, inputs...)
    probs = probabilities(est.discretization, inputs...)
    return association(est.definition, probs)
end

function association(definition::ConditionalEntropyTsallisFuruichi, pxy::Probabilities{T,2}) where {T}
    (; base, q) = definition
    Nx, Ny = size(pxy)
    if q == 1
        return association(ConditionalEntropyShannon(; base), pxy)
    end
    py = marginal(pxy, dims=2)
    ce = 0.0
    qlog = logq0(q)
    for j in 1:Ny
        pyⱼ = py[j]
        for i in 1:Nx
            pxyᵢⱼ = pxy[i, j]
            ce += pxyᵢⱼ^q * qlog(pxyᵢⱼ / pyⱼ)
        end
    end
    ce *= -1.0
    return ce
end


function logq0(q)
    if q == 1.0
        return x -> zero(x)
    else
        return x -> (x^(1 - q) - 1) / (1 - q)
    end
end

