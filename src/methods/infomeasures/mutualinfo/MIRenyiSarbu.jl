export MIRenyiSarbu

"""
    MIRenyiSarbu <: MutualInformation
    MIRenyiSarbu(; base = 2, q = 1.5)

The discrete Rényi mutual information from Sarbu (2014)[^Sarbu2014].

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence.
- Use with [`mutualinfo`](@ref) to compute the raw mutual information. 

## Description

Sarbu (2014) defines discrete Rényi mutual information as the
Rényi ``\\alpha``-divergence between the conditional joint probability mass function
``p(x, y)`` and the product of the conditional marginals, ``p(x) \\cdot p(y)``:

```math
I(X, Y)^R_q =
\\dfrac{1}{q-1}
\\log \\left(
    \\sum_{x \\in X, y \\in Y}
    \\dfrac{p(x, y)^q}{\\left( p(x)\\cdot p(y) \\right)^{q-1}}
\\right)
```

[^Sarbu2014]: Sarbu, S. (2014, May). Rényi information transfer: Partial Rényi transfer
    entropy and partial Rényi mutual information. In 2014 IEEE International Conference
    on Acoustics, Speech and Signal Processing (ICASSP) (pp. 5666-5670). IEEE.

See also: [`mutualinfo`](@ref).
"""
struct MIRenyiSarbu{E <: Renyi} <: MutualInformation{E}
    e::E
    function MIRenyiSarbu(; q = 1.5, base = 2)
        e = Renyi(; q, base)
        new{typeof(e)}(e)
    end
end

function estimate(measure::MIRenyiSarbu, est::Contingency{<:ProbabilitiesEstimator}, x...)
    return estimate(measure, contingency_matrix(est.est, x...))
end

function estimate(measure::MIRenyiSarbu, est::Contingency{<:Nothing}, x...)
    return estimate(measure, contingency_matrix(x...))
end

function estimate(measure::MIRenyiSarbu, pxy::ContingencyMatrix{T, 2}) where {T}
    px = probabilities(pxy, dims = 1)
    py = probabilities(pxy, dims = 2)
    e = measure.e
    q = e.q

    mi = 0.0
    for i in eachindex(px.p)
        for j in eachindex(py.p)
            pxyᵢⱼ = pxy[i, j]
            mi += pxyᵢⱼ^q / ((px[i] * py[j])^(q - 1))
        end
    end
    if mi == 0
        return 0.0
    else
        return (1 / (q - 1) * log(mi)) / log(ℯ, e.base)
    end
end


function estimate(::MIRenyiSarbu, est::ProbabilitiesEstimator, x, y)
    throw(ArgumentError("MIRenyiSarbu not implemented for $(typeof(est))"))
end


function estimate(::MIRenyiSarbu, est::DifferentialEntropyEstimator, args...)
    throw(ArgumentError("MIRenyiSarbu not implemented for $(typeof(est))"))
end
