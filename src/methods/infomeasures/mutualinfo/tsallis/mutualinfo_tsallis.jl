using Entropies: Probabilities, Tsallis
using Entropies: entropy

"""
    Furuichi <: MutualInformationEstimator

A mutual information estimator that computes the Tsallis mutual information as defined in
Furuichi (2006)[^Furuichi2006], based on Tsallis mutual entropy.

[^Furuichi2006]:
    Furuichi, S. Information theoretical properties of Tsallis entropies. J. Math. Phys.
    2006, 47, 023302.
"""
struct Furuichi <: MutualInformationEstimator

end

function mutualinfo(e::Tsallis, est::Furuichi, p::Probabilities, q::Probabilities)
    return entropy(e, p) + entropy(e, q) + joint_entropy(e, p, q)
end
