using ComplexityMeasures: log_with_base

export MIShannon

"""
    MIShannon <: BivariateInformationMeasure
    MIShannon(; base = 2)

The Shannon mutual information ``I_S(X; Y)``.

## Usage

- Use with [`association`](@ref) to compute the raw Shannon mutual information from input data
    using of of the estimators listed below.
- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence using
    the Shannon mutual information.

## Compatible estimators

- [`JointProbabilities`](@ref) (generic)
- [`EntropyDecomposition`](@ref) (generic)
- [`KSG1`](@ref)
- [`KSG2`](@ref)
- [`GaoOhViswanath`](@ref)
- [`GaoKannanOhViswanath`](@ref)
- [`GaussianMI`](@ref)

## Discrete definition

There are many equivalent formulations of discrete Shannon mutual information, meaning that 
it can be estimated in several ways, either using [`JointProbabilities`](@ref)  (double-sum formulation),
[`EntropyDecomposition`](@ref) (three-entropies decomposition), or some dedicated estimator.

### Double sum formulation

Assume we observe samples
``\\bar{\\bf{X}}_{1:N_y} = \\{\\bar{\\bf{X}}_1, \\ldots, \\bar{\\bf{X}}_n \\}`` and
``\\bar{\\bf{Y}}_{1:N_x} = \\{\\bar{\\bf{Y}}_1, \\ldots, \\bar{\\bf{Y}}_n \\}`` from
two discrete random variables ``X`` and ``Y`` with finite supports
``\\mathcal{X} = \\{ x_1, x_2, \\ldots, x_{M_x} \\}`` and
``\\mathcal{Y} = y_1, y_2, \\ldots, x_{M_y}``.
The double-sum estimate is obtained by replacing the double sum

```math
\\hat{I}_{DS}(X; Y) =
 \\sum_{x_i \\in \\mathcal{X}, y_i \\in \\mathcal{Y}} p(x_i, y_j) \\log \\left( \\dfrac{p(x_i, y_i)}{p(x_i)p(y_j)} \\right)
```

where  ``\\hat{p}(x_i) = \\frac{n(x_i)}{N_x}``, ``\\hat{p}(y_i) = \\frac{n(y_j)}{N_y}``,
and ``\\hat{p}(x_i, x_j) = \\frac{n(x_i)}{N}``, and ``N = N_x N_y``.
This definition is used by [`mutualinfo`](@ref) when called with a
[`ContingencyMatrix`](@ref).

### Three-entropies formulation

An equivalent formulation of discrete Shannon mutual information is

```math
I^S(X; Y) = H^S(X) + H_q^S(Y) - H^S(X, Y),
```

where ``H^S(\\cdot)`` and ``H^S(\\cdot, \\cdot)`` are the marginal and joint discrete
Shannon entropies. This definition is used by [`mutualinfo`](@ref) when called with a
[`ProbabilitiesEstimator`](@ref).

## Differential mutual information

One possible formulation of differential Shannon mutual information is

```math
I^S(X; Y) = h^S(X) + h_q^S(Y) - h^S(X, Y),
```

where ``h^S(\\cdot)`` and ``h^S(\\cdot, \\cdot)`` are the marginal and joint
differential Shannon entropies. This definition is used by [`mutualinfo`](@ref) when
called with a [`DifferentialEntropyEstimator`](@ref).

## Estimation

- [Example 1](@ref example_MIShannon_JointProbabilities_ValueBinning): [`JointProbabilities`](@ref) with [`ValueBinning`](@ref) outcome space.
- [Example 2](@ref example_MIShannon_JointProbabilities_UniqueElements): [`JointProbabilities`](@ref) with [`UniqueElements`](@ref) outcome space on string data.
- [Example 3](@ref example_MIShannon_GaussianMI): Dedicated [`GaussianMI`](@ref) estimator.
- [Example 4](@ref example_MIShannon_KSG1): Dedicated [`KSG1`](@ref) estimator.
- [Example 5](@ref example_MIShannon_KSG2): Dedicated [`KSG2`](@ref) estimator.
- [Example 6](@ref example_MIShannon_GaoKannanOhViswanath): Dedicated [`GaoKannanOhViswanath`](@ref) estimator.
- [Example 7](@ref example_MIShannon_EntropyDecomposition_Kraskov): [`EntropyDecomposition`](@ref) with [`Kraskov`](@ref) estimator.
- [Example 8](@ref example_MIShannon_EntropyDecomposition_BubbleSortSwaps): [`EntropyDecomposition`](@ref) with [`BubbleSortSwaps`](@ref).
- [Example 9](@ref example_MIShannon_EntropyDecomposition_Jackknife_ValueBinning): [`EntropyDecomposition`](@ref) with [`Jackknife`](@ref) estimator and [`ValueBinning`](@ref) outcome space.
- [Example 10](@ref example_MIShannon_reproducing_Kraskov): Reproducing Kraskov et al. (2004).
"""
Base.@kwdef struct MIShannon{B} <: MutualInformation
    base::B = 2
end

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:MIShannon}, x, y)
    cts = counts(est.discretization, x, y)
    probs = probabilities(est.discretization, x, y)

    return association(est.definition, probs)
end

function association(definition::MIShannon, pxy::Probabilities{T, 2}) where T
    (; base) = definition
    
    px = marginal(pxy, dims = 1)
    py = marginal(pxy, dims = 2)
    mi = 0.0
    logb = log_with_base(base)
    for i in eachindex(px.p)
        pxᵢ = px[i]
        for j in eachindex(py.p)
            pyⱼ = py[j]
            pxyᵢⱼ = pxy[i, j]
            inner = pxyᵢⱼ / (pxᵢ * pyⱼ)
            if inner != 0.0
                mi += pxyᵢⱼ * logb(inner)
            end
        end
    end
    return mi
end

# ------------------------------------------------
# Mutual information through entropy decomposition
# ------------------------------------------------
function association(est::EntropyDecomposition{<:MIShannon, <:DifferentialInfoEstimator{<:Shannon}}, x, y)
    HX, HY, HXY = marginal_entropies_mi3h_differential(est, x, y)
    mi =  HX + HY - HXY
    return mi
end

function association(est::EntropyDecomposition{<:MIShannon, <:DiscreteInfoEstimator{<:Shannon}, D, P}, x, y) where {D, P}
    HX, HY, HXY = marginal_entropies_mi3h_discrete(est, x, y)
    mi =  HX + HY - HXY
    return mi
end

# ------------------------------------------------
# Pretty printing for decomposition estimators.
# ------------------------------------------------
function decomposition_string(
        definition::MIShannon, 
        est::EntropyDecomposition{<:MIShannon, <:DifferentialInfoEstimator{<:Shannon}}
    )
    return "Iₛ(X, Y) = hₛ(X) + hₛ(Y) - hₛ(X, Y)";
end

function decomposition_string(
        definition::MIShannon, 
        est::EntropyDecomposition{<:MIShannon, <:DiscreteInfoEstimator{<:Shannon}}
    )
    return "Iₛ(X, Y) = Hₛ(X) + Hₛ(Y) - Hₛ(X, Y)";
end
