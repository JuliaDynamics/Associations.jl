export PartialMutualInformation

"""
    PartialMutualInformation <: MultivariateInformationMeasure
    PartialMutualInformation(; base = 2)

The partial mutual information (PMI) measure of conditional association [Zhao2016](@cite).

## Definition

PMI is defined as for variables ``X``, ``Y`` and ``Z`` as

```math
PMI(X; Y | Z) = D(p(x, y, z) || p^{*}(x|z) p^{*}(y|z) p(z)),
```

where ``p(x, y, z)`` is the joint distribution for ``X``, ``Y`` and ``Z``, and
``D(\\cdot, \\cdot)`` is the extended Kullback-Leibler divergence from
``p(x, y, z)`` to ``p^{*}(x|z) p^{*}(y|z) p(z)``. See [Zhao2016](@citet) for details.

## Estimation

The PMI is estimated by first estimating a 3D probability mass function using 
[`probabilities`](@ref), then computing ``PMI(X; Y | Z)`` from those probaiblities.

## Properties

For the discrete case, the following identities hold in theory (when estimating PMI, they
may not).

- `PMI(X, Y, Z) >= CMI(X, Y, Z)` (where CMI is the Shannon CMI). Holds in theory, but
    when estimating PartialMutualInformation, the identity may not hold.
- `PMI(X, Y, Z) >= 0`. Holds both in theory and when estimating using discrete estimators.
- `X ⫫ Y | Z => PMI(X, Y, Z) = CMI(X, Y, Z) = 0` (in theory, but not necessarily for
    estimation).
"""
Base.@kwdef struct PartialMutualInformation <: MultivariateInformationMeasure
    base::Real = 2
end

min_inputs_vars(::PartialMutualInformation) = 3
max_inputs_vars(::PartialMutualInformation) = 3


# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(definition::PartialMutualInformation, pxyz::Probabilities{T, 3}) where T
    dx, dy, dz = size(pxyz)
    px = marginal(pxyz, dims = [1])
    py = marginal(pxyz, dims = [2])
    pz = marginal(pxyz, dims = [3])
    pxz = marginal(pxyz, dims = [1, 3])
    pyz = marginal(pxyz, dims = [2, 3])
    
    # Precompute p⭐x⎸z and p⭐y⎸z sums
    p⭐x⎸z = zeros(dx, dz)
    p⭐y⎸z = zeros(dy, dz)
    
    for i in 1:dx
        for k in 1:dz
            sum_j = 0.0
            for j in 1:dy
                pyⱼ = py[j]
                pyⱼzₖ = pyz[j, k]
                if pyⱼzₖ > 0
                    sum_j += (pxyz[i, j, k] / pyⱼzₖ) * pyⱼ
                end
            end
            p⭐x⎸z[i, k] = sum_j
        end
    end

    for j in 1:dy
        for k in 1:dz
            sum_i = 0.0
            for i in 1:dx
                pxᵢ = px[i]
                pxᵢzₖ = pxz[i, k]
                if pxᵢzₖ > 0
                    sum_i += (pxyz[i, j, k] / pxᵢzₖ) * pxᵢ
                end
            end
            p⭐y⎸z[j, k] = sum_i
        end
    end

    # Compute PMI
    pmi = 0.0
    logb = log_with_base(definition.base)
    for k in 1:dz
        pzₖ = pz[k]
        for j in 1:dy
            p⭐yⱼ⎸zₖ = p⭐y⎸z[j, k]
            if p⭐yⱼ⎸zₖ > 0
                for i in 1:dx
                    p⭐xᵢ⎸zₖ = p⭐x⎸z[i, k]
                    if (p⭐xᵢ⎸zₖ > 0)
                        pxᵢyⱼzₖ = pxyz[i, j, k]
                        if pxᵢyⱼzₖ > 0
                            pxᵢyⱼ⎸zₖ = pxᵢyⱼzₖ / pzₖ
                            logratio = logb(pxᵢyⱼ⎸zₖ  / (p⭐xᵢ⎸zₖ * p⭐yⱼ⎸zₖ))
                            if logratio > 0
                                pmi += pxᵢyⱼzₖ * logratio
                            end
                        end
                    end
                end
            end
        end
    end
    return pmi
end

function association(est::JointProbabilities{PartialMutualInformation}, x, y, z)
    pxyz = probabilities(est.discretization, x, y, z)
    return association(est.definition, pxyz)
end