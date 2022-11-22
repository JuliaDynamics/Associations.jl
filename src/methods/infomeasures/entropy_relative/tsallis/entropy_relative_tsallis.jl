using Entropies: Probabilities
using Entropies: Tsallis

"""
    Tsallis98 <: RelativeEntropyEstimator

A relative entropy estimator that computes the Tsallis relative entropy as defined in
Tsallis (1998)[^Tsallis1998].

[^Tsallis1998]:
    Tsallis, C. Generalized entropy-based criterion for consistent testing.
    Phys. Rev. E 1998, 58, 479-487.
"""
struct Tsallis98 <: RelativeEntropyEstimator end

function entropy_relative(e::Tsallis, est::Tsallis98, p::Probabilities, q::Probabilities)
    q = e.q
    return 1 / (q - 1) * (1 - sum(pᵢ^q / qᵢ^(q - 1) for (pᵢ, qᵢ) in zip(p, q)))
end
