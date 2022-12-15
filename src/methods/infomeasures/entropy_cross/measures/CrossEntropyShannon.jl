"""
    CrossEntropyShannon <: CrossEntropy
    CrossEntropyShannon(; base = 2)

`CrossEntropyShannon` is a directive to compute the discrete Shannon cross-entropy
(or divergence) to the given `base` between random variables
``X \\in \\mathbb{R}^{d_X}`` and ``Y \\in \\mathbb{R}^{d_Y}`` using the formula.

## Supported definitions

See also: [`entropy_cross`](@ref).
"""
struct CrossEntropyShannon{E <: Renyi} <: CrossEntropy
    e::E
    function CrossEntropyShannon(; base = 2)
            e = Shannon(; base)
        @assert e.q ≈ 1 || error("CrossEntropyShannon not defined for q = $(e.q)")
        new{typeof(e)}(e)
    end
end
