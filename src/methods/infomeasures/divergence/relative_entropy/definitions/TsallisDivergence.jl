export TsallisDivergence
export TsallisDivergenceDifferential

"""
    TsallisDivergence <: DivergenceDefinition
    TsallisDivergence()

`TsallisDivergence` gives the original definition of Tsallis relative entropy
from Tsallis (1998)[^Tsallis1998].

## Description

For a discrete sample space ``\\Omega`` and probability mass functions
``a(x) : \\Omega \\to [0, 1]`` and ``b(x) : \\Omega \\to [0, 1]``,
the Tsallis relative entropy (divergence) is given by

```math
D_q^T(A || B) = - \\sum_{i = 1}^n a_i^q \\log_q \\dfrac{b_i}{a_i}.
```

For `q = 1`, the Shannon relative entropy (KL divergence) is obtained, and ``log_q`` is
the q-logarithm defined in Tsallis (1998).

[^Tsallis1998]:
    Tsallis, C. (1998). Generalized entropy-based criterion for consistent testing.
    Physical Review E, 58(2), 1442.
"""
struct TsallisDivergence <: DivergenceDefinition end

function divergence(def::TsallisDivergence, p::Probabilities, q::Probabilities)
    q = e.q
    return 1 / (q - 1) * (1 - sum(pᵢ^q / qᵢ^(q - 1) for (pᵢ, qᵢ) in zip(p, q)))
end
