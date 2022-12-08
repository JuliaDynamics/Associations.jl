using Entropies: Probabilities, Tsallis
using Entropies: entropy

"""
    TsallisMIFuruichi <: MutualInformationEstimator
    TsallisMIFuruichi(est::ProbabilitiesEstimator)

`TsallisMIFuruichi` is an estimator of the Tsallis mutual "information",
which they Furuichi (2006)[^Furuichi2006] denotes the Tsallis *mutual entropy*
(see Furuichi (2006) for reasoning on naming the method).

[^Furuichi2006]:
    Furuichi, S. Information theoretical properties of Tsallis entropies. J. Math. Phys.
    2006, 47, 023302.
"""
struct TsallisMIFuruichi{P} <: MutualInformationEstimator
    est::P
end

function mutualinfo(e::Tsallis, est::TsallisMIFuruichi, x, y)
    return entropy(e, x) + entropy(e, y) + entropy(e, Dataset(x, y))
end
