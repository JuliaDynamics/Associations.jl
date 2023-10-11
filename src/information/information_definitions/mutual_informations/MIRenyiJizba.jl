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
I_q^{R_{J}}(X; Y) = S_q^{R}(X) + S_q^{R}(Y) - S_q^{R}(X, Y),
```

where ``S_q^{R}(\\cdot)`` and ``S_q^{R}(\\cdot, \\cdot)`` the [`Rényi`](@ref) entropy and
the joint Rényi entropy.
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


function information(est::EntropyDecomposition{<:MIRenyiJizba}, x, y)
    throw(ArgumentError("MIRenyiSarbu not implemented for $(typeof(est).name.name)"))
end
