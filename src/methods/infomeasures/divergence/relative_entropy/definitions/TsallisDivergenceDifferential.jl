
export TsallisDivergenceDifferential

"""
    TsallisDivergenceDifferential <: DivergenceDefinition
    TsallisDivergenceDifferential()

`TsallisDivergenceDifferential` is a directive to compute the Tsallis differential
relative entropy (divergence) defined in Póczos & Schneider (2011)[^Póczos2011].

## Description

Let ``\\mathbb{P}`` and ``\\mathbb{Q}`` be continuous probability measures
with density functions ``p(x)`` and ``q(x)`` (``x \\in M_0 \\subseteq \\mathcal{R}^D``),
i.e. ``p: M_0 \\subseteq \\mathcal{R}^D \\mapsto \\mathbb{R}`` and
``q: M_0 \\subseteq \\mathcal{R}^D \\mapsto \\mathbb{R}``,
with respect to the Lebesque measure ``\\mu``, and ``dx := \\mu(dx)``.

Let ``q\\in \\mathbb{R}`` with ``q \\neq 1``. Assuming the following
integral exists:

```math
V_{q}(\\mathbb{P} || \\mathbb{Q}) = \\int_{M_0} p^{q}(x)q^{1 - q}(x) dx.
```

Then the differential Tsallis divergence is given by

```math
T_{q}(\\mathbb{P} || \\mathbb{Q}) =
\\dfrac{1}{q - 1}
\\left(
    V_{q}(\\mathbb{P} || \\mathbb{Q}) - 1
\\right).
```

[^Póczos2011]:
    Póczos, B., & Schneider, J. (2011, June). On the estimation of alpha-divergences. In
    Proceedings of the Fourteenth International Conference on Artificial Intelligence and
    Statistics (pp. 609-617). JMLR Workshop and Conference Proceedings.
"""
struct TsallisDivergenceDifferential <: DivergenceDefinition end
