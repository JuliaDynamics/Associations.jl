export PMI
export pmi

"""
    PMI <: AssociationMeasure
    PMI(; base = 2)

The partial mutual information (PMI) measure of association (Zhao et al., 2016)[^Zhao2016].

## Definition

PMI is defined as for variables ``X``, ``Y`` and ``Z`` as

```math
PMI(X; Y | Z) = D(p(x, y, z) || p^{*}(x|z) p^{*}(y|z) p(z)),
```

where ``p(x, y, z)`` is the joint distribution for ``X``, ``Y`` and ``Z``, and
``D(\\cdot, \\cdot)`` is the extended Kullback-Leibler divergence from
``p(x, y, z)`` to ``p^{*}(x|z) p^{*}(y|z) p(z)``. See Zhao et al. (2016) for details.

## Estimation

PMI can be estimated using any [`ProbabilitiesEstimator`](@ref) that implements
[`marginal_encodings`](@ref). This allows estimation of 3D contingency matrices, from
which relevant probabilities for the PMI formula are extracted. See also [`pmi`](@ref).

## Properties

For the discrete case, the following identities hold in theory (when estimating PMI, they
may not).

- `PMI(X, Y, Z) >= CMI(X, Y, Z)` (where CMI is the Shannon CMI). Holds in theory, but
    when estimating PMI, the identity may not hold.
- `PMI(X, Y, Z) >= 0`. Holds both in theory and for estimation using
    [`ProbabilitiesEstimator`](@ref)s.
- `X ⫫ Y | Z => PMI(X, Y, Z) = CMI(X, Y, Z) = 0` (in theory, but not necessarily for
    estimation).

[^Zhao2016]:
    Zhao, J., Zhou, Y., Zhang, X., & Chen, L. (2016). Part mutual information for
    quantifying direct associations in networks. Proceedings of the National Academy
    of Sciences, 113(18), 5130-5135.
"""
Base.@kwdef struct PMI <: AssociationMeasure
    base::Real = 2
end

"""
    pmi([measure::CMI], est::ProbabilitiesEstimator, x, y, z) → pmi_est::Real ∈ [0, a)

Estimate the part mutual information ([`PMI`](@ref); Zhao et al. (2016)[^Zhao2016]).

If `measure` is not given, then the default is `PMI(; base = 2)`.
With a [`ProbabilitiesEstimator`](@ref), the returned `pmi_est` is guaranteed to be
non-negative.

## Estimators

| Estimator                     | Principle           | [`PMI`](@ref) |
| ----------------------------- | ------------------- | :-----------: |
| [`CountOccurrences`](@ref)    | Frequencies         |          ✓    |
| [`ValueHistogram`](@ref)      | Binning (histogram) |          ✓    |
| [`OrdinalPatterns`](@ref) | Ordinal patterns    |          ✓    |
| [`Dispersion`](@ref)          | Dispersion patterns |          ✓    |

[^Zhao2016]:
    Zhao, J., Zhou, Y., Zhang, X., & Chen, L. (2016). Part mutual information for
    quantifying direct associations in networks. Proceedings of the National Academy
    of Sciences, 113(18), 5130-5135.
"""
function pmi(measure::PMI, x...)
    return estimate(measure, x...)
end

function pmi(x...)
    return estimate(PMI(), x...)
end

function estimate(measure::PMI, est::Contingency{<:ProbabilitiesEstimator}, x...)
    return estimate(measure, contingency_matrix(est.est, x...))
end

function estimate(measure::PMI, est::Contingency{<:Nothing}, x...)
    return estimate(measure, contingency_matrix(CountOccurrences(), x...))
end

# We explicitly need to construct a contingency matrix, because unlike for e.g. CMI,
# there's no obvious way to rewrite PMI in terms of sums of entropies.
function estimate(measure::PMI, est::ProbabilitiesEstimator, x...)
    return estimate(measure, Contingency(est), x...)
end

function estimate(
        measure::PMI,
        pxyz::ContingencyMatrix{T, 3}) where T

    # The sums go over *states*, so these are what we iterate over.
    dx, dy, dz = size(pxyz)
    px = probabilities(pxyz, dims = [1])
    py = probabilities(pxyz, dims = [2])
    pz = probabilities(pxyz, dims = [3])
    pyz = probabilities(pxyz, dims = [2, 3])
    pxz = probabilities(pxyz, dims = [1, 3])

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
