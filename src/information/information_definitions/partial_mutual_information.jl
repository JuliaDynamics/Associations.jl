export PartialMutualInformation

"""
    PartialMutualInformation <: MultivariateInformationMeasure
    PartialMutualInformation(; base = 2)

The partial mutual information (PMI) measure of conditional association [Zhao2016](@cite).

## Definition

PMI is defined as for variables ``X``, ``Y`` and ``Z`` as

```math
PartialMutualInformation(X; Y | Z) = D(p(x, y, z) || p^{*}(x|z) p^{*}(y|z) p(z)),
```

where ``p(x, y, z)`` is the joint distribution for ``X``, ``Y`` and ``Z``, and
``D(\\cdot, \\cdot)`` is the extended Kullback-Leibler divergence from
``p(x, y, z)`` to ``p^{*}(x|z) p^{*}(y|z) p(z)``. See Zhao et al. (2016) for details.

## Estimation

The PMI is estimated by first estimating a 3D probability mass function using 
[`probabilities`](@ref), then computing ``PartialMutualInformation(X; Y | Z)`` from those probaiblities.

## Properties

For the discrete case, the following identities hold in theory (when estimating PMI, they
may not).

- `PartialMutualInformation(X, Y, Z) >= CMI(X, Y, Z)` (where CMI is the Shannon CMI). Holds in theory, but
    when estimating PartialMutualInformation, the identity may not hold.
- `PartialMutualInformation(X, Y, Z) >= 0`. Holds both in theory and when estimating using discrete estimators.
- `X ⫫ Y | Z => PartialMutualInformation(X, Y, Z) = CMI(X, Y, Z) = 0` (in theory, but not necessarily for
    estimation).
"""
Base.@kwdef struct PartialMutualInformation <: MultivariateInformationMeasure
    base::Real = 2
end

min_inputs_vars(::PartialMutualInformation) = 3
max_inputs_vars(::PartialMutualInformation) = 3

function information(definition::PartialMutualInformation, pxyz::Probabilities{T, 3}) where T
    dx, dy, dz = size(pxyz)
    px = marginal(pxyz, dims = [1])
    py = marginal(pxyz, dims = [2])
    pz = marginal(pxyz, dims = [3])
    pyz = marginal(pxyz, dims = [2, 3])
    pxz = marginal(pxyz, dims = [1, 3])

    pmi = 0.0
    logb = log_with_base(measure.base)
    for i in 1:dx
        for j in 1:dy
            for k in 1:dz
                pzₖ = pz[k]
                pyzⱼₖ = pyz[j, k]
                pxzᵢₖ = pxz[i, k]
                pxyzᵢⱼₖ = pxyz[i, j, k]
                sy = sum(pyzⱼₖ > 0 ? py[j] * (pxyz[i, j, k]/ pyz[j, k]) : 0 for j = 1:dy)
                sx = sum(pxzᵢₖ > 0 ? px[i] * (pxyz[i, j, k] / pxz[i, k]) : 0 for i = 1:dx)

                sxy = sy * sx
                pxy_z = pxyzᵢⱼₖ / pzₖ
                if sxy > 0 && pxy_z > 0
                    pmi += pxyzᵢⱼₖ * logb(pxy_z / (sy * sx))
                end
            end
        end
    end
    return pmi
end