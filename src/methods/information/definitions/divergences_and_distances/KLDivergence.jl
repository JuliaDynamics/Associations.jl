export KLDivergence

"""
    KLDivergence <: DivergenceOrDistance

The Kullback-Leibler (KL) divergence.

## Usage 

- Use with [`association`](@ref) to compute the compute the KL-divergence between two 
    pre-computed probability distributions, or from raw data using one of the estimators
    listed below.

## Compatible estimators

- [`JointDistanceDistribution`](@ref)

## Estimators 

- [`JointProbabilities`](@ref).

## Description

The KL-divergence between two probability distributions
``P_X = (p_x(\\omega_1), \\ldots, p_x(\\omega_n))`` and
``P_Y = (p_y(\\omega_1), \\ldots, p_y(\\omega_m))``, both defined over the same
[`OutcomeSpace`](@ref) ``\\Omega = \\{\\omega_1, \\ldots, \\omega_n \\}``, is defined as

```math
D_{KL}(P_Y(\\Omega) || P_Y(\\Omega)) =
\\sum_{\\omega \\in \\Omega} p_x(\\omega) \\log\\dfrac{p_x(\\omega)}{p_y(\\omega)}
```

## Implements

- [`association`](@ref). Used to compute the KL-divergence between two pre-computed
    probability distributions. If used with [`RelativeAmount`](@ref), the KL divergence may
    be undefined to due some outcomes having zero counts. Use some other
    [`ProbabilitiesEstimator`](@ref) like [`BayesianRegularization`](@ref) to ensure
    all estimated probabilities are nonzero.

!!! note 
    Distances.jl also defines `KLDivergence`. Quality it if you're loading both 
    packages, i.e. do `association(CausalityTools.KLDivergence(), x, y)`.

## Estimation

- [Example 1](@ref example_KLDivergence_precomputed_probabilities): From precomputed probabilities
- [Example 2](@ref example_KLDivergence_JointProbabilities_OrdinalPatterns): 
    [`JointProbabilities`](@ref) with [`OrdinalPatterns`](@ref) outcome space
"""
struct KLDivergence{B} <: DivergenceOrDistance
    base::B
end
KLDivergence(; base = 2) = KLDivergence(base)

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:KLDivergence}, x, y)
    # Dispatch to generic method in `divergences_and_distances.jl` with 2D `Probabilities`
    probs = probabilities(est.discretization, x, y)
    return association(est.definition, probs)
end

function association(measure::KLDivergence, px::Probabilities, py::Probabilities)
    size_match(measure, px, py)
    return sum(pxᵢ * log(measure.base, pxᵢ / pyᵢ) for (pxᵢ, pyᵢ) in zip(px, py))
end
