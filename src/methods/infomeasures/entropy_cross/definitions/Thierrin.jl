export Thierrin

"""
    Thierrin <: CrossEntropyDefinition
    Thierrin()

`Thierrin` is a directive to compute Rényi cross-entropy as defined by
Thierrin et al. (2022):

```math
C_{q}^R(\\mathcal{P}; \\mathcal{Q}) :=
\\dfrac{1}{1-q}
\\log \\sum_{x \\in S} {\\left( p(x)q(x)^{q-1} \\right)}.
```

[^Thierrin2022]:
    Thierrin, F. C., Alajaji, F., & Linder, T. (2022). Rényi Cross-Entropy Measures for
    Common Distributions and Processes with Memory. Entropy, 24(10), 1417.
"""
struct Thierrin <: CrossEntropyDefinition end

function estimate(def::Thierrin, measure::CrossEntropyRenyi,
        est::ProbabilitiesEstimator, x, y)
    e = measure.e
    q, base = e.q, e.base
    px = probabilities(est, x)
    py = probabilities(est, y)
    if q ≈ 1
        return -sum(pxᵢ * log0(base, pyᵢ) for (pxᵢ, pyᵢ) in zip(px, py))
    else
        return -sum(pxᵢ * log0(base, pyᵢ^(q-1)) for (pxᵢ, pyᵢ) in zip(px, py))
    end
end

estimate(measure::CrossEntropyRenyi, est, x, y) = estimate(Thierrin(), measure, est, x, y)
