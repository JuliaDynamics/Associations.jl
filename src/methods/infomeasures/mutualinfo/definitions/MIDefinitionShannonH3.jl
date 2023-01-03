export MIDefinitionShannonH3

"""
    MIDefinitionShannonH3 <: MutualInformationDefinition

The 3-entropies formulation of Shannon mutual information ``I^S(X; Y)``

This is the maximum likelihood estimate for discrete mutual information (Paninski, 2003).

## Description

Assume we observe samples
``\\bar{\\bf{X}}_{1:n} = \\{\\bar{\\bf{X}}_1, \\ldots, \\bar{\\bf{X}}_n \\}`` and
``\\bar{\\bf{Y}}_{1:n} = \\{\\bar{\\bf{Y}}_1, \\ldots, \\bar{\\bf{Y}}_n \\}`` from
two discrete random variables ``X`` and ``Y`` with finite supports ``\\mathcal{X}`` and
``\\mathcal{Y}``. The H3-definition of mutual information is

- Continuous case: ``I^S(X; Y) = H^S(X) + H_q^S(Y) - H^S(X, Y)``
- Discrete case: ``I^S(X; Y) = H^S(X) + H_q^S(Y) - H^S(X, Y)``,

Here, ``H^S(\\cdot)`` and ``H^S(\\cdot, \\cdot)`` are the marginal and joint discrete
Shannon entropies, and ``h^S(\\cdot)`` and ``h^S(\\cdot, \\cdot)`` are the corresponding
differential entropies. To compute ``I^S(X; Y)``, we simply compute each marginal
entropy separately and sum them according to the formuals above.

## Compatibility (continuous estimators)

This definition is compatible with all [`DifferentialEntropyEstimator`](@ref) that
accept multivariate inputs.

- [`Kraskov`](@ref)
- [`KozachenkoLeonenko`](@ref)
- [`Gao`](@ref)
- [`Goria`](@ref)
- [`Zhu`](@ref)
- [`ZhuSingh`](@ref)
- [`LeonenkoProzantoSavani`](@ref)
- [`Lord`](@ref)

## Compatibility (discrete estimators)

This definition is compatible with any [`ProbabilitiesEstimator`](@ref) that accepts
multivariate input data (required for estimation of the joint probabilities).
In addition, we provide specialized methods for the following estimators:

- **[`SymbolicPermutation`](@ref)**. `X` and `Y` are encoded into their ordinal patterns
     separately, and the joint space `XY` is simply the concatenation of the
     ordinal pattern vectors.
- **[`Dispersion`](@ref)**. `X` and `Y` are encoded into their dispersion patterns
    separately, and the joint space `XY` is simply the concatenation of the
    ordinal pattern vectors.
- **[`ValueHistogram`](@ref)**. Bin visitation frequencies are counted in the joint space
    `XY`, then marginal probabilities are obtained from the joint bin visits.
"""
struct MIDefinitionShannonH3 <: MIH3

end
