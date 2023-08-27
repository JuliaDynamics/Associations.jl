export CEShannon

"""
    CEShannon <: ConditionalEntropy
    CEShannon(; base = 2,)

The[`Shannon`](@ref) conditional entropy measure.

## Discrete definition

### Sum formulation

The conditional entropy between discrete random variables
``X`` and ``Y`` with finite ranges ``\\mathcal{X}`` and ``\\mathcal{Y}`` is defined as

```math
H^{S}(X | Y) = -\\sum_{x \\in \\mathcal{X}, y \\in \\mathcal{Y}} = p(x, y) \\log(p(x | y)).
```

This is the definition used when calling [`entropy_conditional`](@ref) with a
[`ContingencyMatrix`](@ref).

### Two-entropies formulation

Equivalently, the following difference of entropies hold

```math
H^S(X | Y) = H^S(X, Y) - H^S(Y),
```

where ``H^S(\\cdot`` and ``H^S(\\cdot | \\cdot)`` are the [`Shannon`](@ref) entropy and
Shannon joint entropy, respectively. This is the definition used when calling
[`entropy_conditional`](@ref) with a [`OutcomeSpace`](@ref).

## Differential definition

The differential conditional Shannon entropy is analogously defined as

```math
H^S(X | Y) = h^S(X, Y) - h^S(Y),
```

where ``h^S(\\cdot`` and ``h^S(\\cdot | \\cdot)`` are the [`Shannon`](@ref)
differential entropy and Shannon joint differential entropy, respectively. This is the
definition used when calling [`entropy_conditional`](@ref) with a
[`DifferentialInfoEstimator`](@ref).
"""
struct CEShannon{E} <: ConditionalEntropy
    e::E
    function CEShannon(; base = 2)
        e = MLEntropy(Shannon(; base))
        new{typeof(e)}(e)
    end
end

function estimate(measure::CEShannon, c::ContingencyMatrix{T, 2}) where {T}
    e = measure.e.definition
    Nx, Ny = size(c)

    py = probabilities(c, dims = 2)
    pxy = probabilities(c)
    ce = 0.0
    log0 = log_with_base(e.base)
    for j in 1:Ny
        pyⱼ = py[j]
        for i in 1:Nx
            pxyᵢⱼ = pxy[i, j]
            if pxyᵢⱼ != 0.0
                ce += pxyᵢⱼ * log0(pxyᵢⱼ / pyⱼ)
            end
        end
    end
    return -ce
end

function estimate(measure::CEShannon, est::DiscreteOrDifferentialInfoEstimator, x, y)
    HY, HXY = marginal_entropies_ce2h(measure, est, x, y)
    return HXY - HY
end
