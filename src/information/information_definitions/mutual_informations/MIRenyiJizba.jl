using ComplexityMeasures: Renyi

export MIRenyiJizba

"""
    MIRenyiJizba <: <: BivariateInformationMeasure
    MIRenyiJizba(; q = 1.5, base = 2)

The Rényi mutual information ``I_q^{R_{J}}(X; Y)`` defined in [Jizba2012](@cite).

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence.
- Use with [`mutualinfo`](@ref) to compute the raw mutual information.

## Definition

```math
I_q^{R_{J}}(X; Y) = H_q^{R}(X) + H_q^{R}(Y) - H_q^{R}(X, Y),
```

where ``H_q^{R}(\\cdot)`` is the [`Rényi`](@ref) entropy.

## 

## Examples

For categorical data, we can use the `JointProbabilities` estimator.

```julia
using CausalityTools, Random; rng = MersenneTwister(1234)
x = rand(rng, ["a", "b", "c"], 200);
y = rand(rng, ["hello", "yoyo", "heyhey"], 200);
est = JointProbabilities(MIRenyiJizba(), UniqueElements())
information(est, x, y)
```

For numeric data, we can use the `EntropyDecomposition` estimator.

```julia
x = randn(rng, 50);
y = randn(rng, 50);
def = MIRenyiJizba()
est_diff = EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi(), k=3))
information(est_diff, x, y) 

est_disc = EntropyDecomposition(def, PlugIn(Renyi()), ValueBinning(2));
information(est_disc, x, y)
````
"""
Base.@kwdef struct MIRenyiJizba{B, Q} <: MutualInformation
    base::B = 2
    q::Q = 1.5
end

function information(definition::MIRenyiJizba, pxy::Probabilities{T, 2}) where T
    (; base, q) = definition

    px = marginal(pxy, dims = 1)
    py = marginal(pxy, dims = 2)
    
    logb = log_with_base(base)
    num = 0.0
    den = 0.0
    for i in eachindex(px.p)
        for j in eachindex(py.p)
            num += px[i]^q * py[j]^q
            den += pxy[i, j]^q
        end
    end
    if den != 0
        mi = logb(num / den)
    else
        mi = 0.0
    end

    return (1 / (1 / q)) * mi
end

# --------------------------------------------------------------
# `MIRenyiJizba` through entropy decomposition.
# Eq. 24 in
# Jizba, P., Lavička, H., & Tabachová, Z. (2021). Rényi Transfer Entropy Estimators for
# Financial Time Series. Engineering Proceedings, 5(1), 33.
# --------------------------------------------------------------
function information(est::EntropyDecomposition{<:MIRenyiJizba, <:DifferentialInfoEstimator{<:Renyi}}, x, y)
    HX, HY, HXY = marginal_entropies_mi3h_differential(est, x, y)
    mi =  HX + HY - HXY
    return mi
end

function information(est::EntropyDecomposition{<:MIRenyiJizba, <:DiscreteInfoEstimator{<:Renyi}}, x, y)
    HX, HY, HXY = marginal_entropies_mi3h_discrete(est, x, y)
    mi =  HX + HY - HXY
    return mi
end

# ------------------------------------------------
# Pretty printing for decomposition estimators.
# ------------------------------------------------
function decomposition_string(
        definition::MIRenyiJizba, 
        est::EntropyDecomposition{<:MIRenyiJizba, <:DifferentialInfoEstimator{<:Renyi}}
    )
    return "Iᵣⱼ(X, Y) = hᵣ(X) + hᵣ(Y) - hᵣ(X, Y)"
end

function decomposition_string(
        definition::MIRenyiJizba, 
        est::EntropyDecomposition{<:MIRenyiJizba, <:DiscreteInfoEstimator{<:Renyi}}
    )
    return "Iᵣⱼ(X, Y) = Hᵣ(X) + Hᵣ(Y) - Hᵣ(X, Y)"
end
