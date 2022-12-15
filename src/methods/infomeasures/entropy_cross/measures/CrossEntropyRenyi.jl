"""
    CrossEntropyRenyi <: CrossEntropy
    CrossEntropyRenyi(; base = 2, q = 1.5)

`CrossEntropyRenyi` is a directive to compute the discrete Rényi cross-entropy
(or divergence) to the given `base` between random variables
``X \\in \\mathbb{R}^{d_X}`` and ``Y \\in \\mathbb{R}^{d_Y}`` using the formula.

## Supported definitions

See also: [`entropy_cross`](@ref).
"""
struct CrossEntropyRenyi{E <: Renyi} <: CrossEntropy
    e::E
    function CrossEntropyRenyi(; base = 2, q = 1.5)
            e = Renyi(; base, q)
        new{typeof(e)}(e)
    end
end
