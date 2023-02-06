export TERenyiJizba

"""
    TERenyiJizba() <: TransferEntropy

The Rényi transfer entropy from Jizba et al. (2012)[^Jizba2012].

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise
    and conditional dependence.
- Use with [`transferentropy`](@ref) to compute the raw transfer entropy.

## Description

The transfer entropy from source ``S`` to target ``T``, potentially
conditioned on ``C`` is defined as

```math
\\begin{align*}
TE(S \\to T) &:= I_q^{R_J}(T^+; S^- | T^-) \\\\
TE(S \\to T | C) &:= I_q^{R_J}(T^+; S^- | T^-, C^-),
\\end{align*},
```
where ``I_q^{R_J}(T^+; S^- | T^-)`` is Jizba et al. (2012)'s definition of
conditional mutual information ([`CMIRenyiJizba`](@ref)).
The variables ``T^+``, ``T^-``,
``S^-`` and ``C^-`` are described in the docstring for [`transferentropy`](@ref).

## Compatible estimators

- **[`ProbabilitiesEstimator`](@ref)**: Any probabilities estimator that accepts
    multivariate input data or has an implementation for [`marginal_encodings`](@ref).
    Transfer entropy is computed a sum of marginal (discrete) entropy estimates.
    Example: [`ValueHistogram`](@ref).
- **[`DifferentialEntropyEstimator`](@ref)**. Any differential entropy
    estimator that accepts multivariate input data.
    Transfer entropy is computed a sum of marginal differential entropy estimates.
    Example: [`Kraskov`](@ref).

[^Jizba2012]:
    Jizba, P., Kleinert, H., & Shefaat, M. (2012). Rényi’s information transfer between
    financial time series. Physica A: Statistical Mechanics and its Applications, 391(10),
    2971-2989.
"""
struct TERenyiJizba{E <: Renyi, EMB} <: TransferEntropy{E, EMB}
    e::E
    embedding::EMB
    function TERenyiJizba(; base = 2, q = 1.5, embedding::EMB = EmbeddingTE()) where EMB
        e = Renyi(; base = base, q = q)
        return new{typeof(e), EMB}(e, embedding)
    end
    function TERenyiJizba(e::E; embedding::EMB = EmbeddingTE()) where {E <: Renyi, EMB}
        return new{E, EMB}(e, embedding)
    end
end

"""
    escort_distribution(probs, i::Int, q::Real)

The escort distribution for a probability distribution `probs`. For `q > 1`, the
escort distribution emphasises more probable events and de-emphasises more improbable
events. For `q < 1`, the situation is reversed.

```math
\\text{esc}_q(x) = \\dfrac{p^q(x)}{\\sum_{x \\in \\mathcal{X}} p^q(x)}
```
"""
function escort_distribution(probs, i::Int, q)
    return probs[i]^q / sum(probs .^ q)
end
