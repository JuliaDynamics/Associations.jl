export RenyiDivergence

"""
    RenyiDivergence <: DivergenceOrDistance
    RenyiDivergence(q; base = 2)

The Rényi divergence of positive order `q`.

## Usage 

- Use with [`association`](@ref) to compute the compute the Rényi divergence between two 
    pre-computed probability distributions, or from raw data using one of the estimators
    listed below.

## Compatible estimators

- [`JointDistanceDistribution`](@ref)

## Description

The Rényi divergence between two probability distributions
``P_X = (p_x(\\omega_1), \\ldots, p_x(\\omega_n))`` and
``P_Y = (p_y(\\omega_1), \\ldots, p_y(\\omega_m))``, both defined over the same
[`OutcomeSpace`](@extref ComplexityMeasures.OutcomeSpace) ``\\Omega = \\{\\omega_1, \\ldots, \\omega_n \\}``, is defined as
[vanErven2014](@citet).

```math
D_{q}(P_Y(\\Omega) || P_Y(\\Omega)) =
\\dfrac{1}{q - 1} \\log \\sum_{\\omega \\in \\Omega}p_x(\\omega)^{q}p_y(\\omega)^{1-\\alpha}
```

## Implements

- [`information`](@ref). Used to compute the Rényi divergence between two pre-computed
    probability distributions. If used with [`RelativeAmount`](@extref ComplexityMeasures.RelativeAmount), the KL divergence may
    be undefined to due some outcomes having zero counts. Use some other
    [`ProbabilitiesEstimator`](@extref ComplexityMeasures.ProbabilitiesEstimator) like [`BayesianRegularization`](@extref ComplexityMeasures.BayesianRegularization) to ensure
    all estimated probabilities are nonzero.

!!! note 
    Distances.jl also defines `RenyiDivergence`. Quality it if you're loading both 
    packages, i.e. do `association(Associations.RenyiDivergence(), x, y)`.


## Estimation

- [Example 1](@ref example_RenyiDivergence_precomputed_probabilities): From precomputed probabilities
- [Example 2](@ref example_RenyiDivergence_JointProbabilities_OrdinalPatterns): 
    [`JointProbabilities`](@ref) with [`OrdinalPatterns`](@extref ComplexityMeasures.OrdinalPatterns) outcome space
"""
struct RenyiDivergence{Q,B} <: DivergenceOrDistance
    q::Q
    base::B
    function RenyiDivergence(q::Q, base::B) where {Q,B}
        q > 0 || throw(ArgumentError("`q` must be positive. Got $q"))
        new{Q,B}(q, base)
    end
end
RenyiDivergence(; q=0.5, base=2) = RenyiDivergence(q, base)

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:RenyiDivergence}, x, y)
    # Dispatch to generic method in `divergences_and_distances.jl` with 2D `Probabilities`
    probs = probabilities(est.discretization, x, y)
    return association(est.definition, probs)
end

function association(definition::RenyiDivergence, px::Probabilities, py::Probabilities)
    (; base, q) = definition

    if q == Inf
        return maximum(pxᵢ / pyᵢ for (pxᵢ, pyᵢ) in zip(px, py))
    end
    s = 0.0
    for (pxᵢ, pyᵢ) in zip(px, py)
        if pxᵢ != 0.0 && pyᵢ != 0.0
            s += pxᵢ^q * pyᵢ^(1 - q)
        end
    end
    return 1 / (q - 1) * log(base, s)
end
