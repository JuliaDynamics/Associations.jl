using Distributions: Exponential, Normal, MvNormal
using LinearAlgebra: inv, det, transpose

export RenyiDifferentialCrossEntropyThierrin2022
export Exponential, Normal, MvNormal
export entropy_cross

"""
    RenyiDifferentialCrossEntropyThierrin2022 <: CrossEntropyDefinition
    RenyiDifferentialCrossEntropyThierrin2022()

Thierrin et al. (2022) defines discrete Rényi cross entropy as

```math
C_{q}^R(\\mathcal{P}; \\mathcal{Q}) :=
\\dfrac{1}{1-q}
\\log \\sum_{x \\in S} {\\left( p(x)q(x)^{q-1} \\right)}.
```

[^Thierrin2022]:
    Thierrin, F. C., Alajaji, F., & Linder, T. (2022). Rényi Cross-Entropy Measures for
    Common Distributions and Processes with Memory. Entropy, 24(10), 1417.
"""
struct RenyiDifferentialCrossEntropyThierrin2022 <: CrossEntropyDefinition end

function entropy_cross(e::Entropy,
        est::DiscreteCrossEntropy{P, RenyiDifferentialCrossEntropyThierrin2022}, x, y) where P
    px = probabilities(est.est, x)
    py = probabilities(est.est, y)

    q = e.q
    if q ≈ 1
        return -sum(pxᵢ * log0(e.base, pyᵢ) for (pxᵢ, pyᵢ) in zip(px, py))
    else
        return -sum(pxᵢ * log0(e.base, pyᵢ^(q-1)) for (pxᵢ, pyᵢ) in zip(px, py))
    end
end
