export RenyiDivergence

"""
    RenyiDivergence  <: DivergenceDefinition
    RenyiDivergence()

An instruction to compute the Rényi divergence (relative entropy) according to the
original definition in Rényi (1961)[^Rényi1961]'s seminal paper, here stated in terms of
notation from Van Erven et al. (2014)[[^VanErven2014]].

## Description

For a discrete sample space ``\\Omega`` and probability mass functions
``p(x) : \\Omega \\to [0, 1]`` and ``q(x) : \\Omega \\to [0, 1]``,
the Rényi relative entropy (divergence) is, for `q != 1`, given by

```math
D_q(P || Q) = \\dfrac{1}{q - 1} \\log \\sum_{i = 1}^n p_i^q q_i^{1-\\alpha}.
```

[^Rényi1961]:
    Rényi, A. (1961, June). On measures of entropy and information. In Proceedings of the
    fourth Berkeley symposium on mathematical statistics and probability (Vol. 1, No.
    547-561).

[^VanErven2014]:
    Van Erven, T., & Harremos, P. (2014). Rényi divergence and Kullback-Leibler divergence.
    IEEE Transactions on Information Theory, 60(7), 3797-3820.
"""
struct RenyiDivergence <: DivergenceDefinition end
