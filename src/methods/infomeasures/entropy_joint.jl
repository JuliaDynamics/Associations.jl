using ComplexityMeasures: EntropyDefinition

export entropy_joint
export JointEntropy
export JointEntropyRenyi
export JointEntropyShannon
export JointEntropyTsallis

"""
    entropy_joint(e::EntropyDefinition, x, y)
    entropy_joint(e::EntropyDefinition, c::ContingencyMatrix)

Compute the joint entropy of type `e` (e.g. [`Shannon`](@ref)) of the input variables
`x` and `y`, or from the pre-computed contingency matrix `c` (see
[`ContingencyMatrix`](@ref)).

## Discrete definitions

Given two two discrete random variables ``X`` and ``Y`` with ranges ``\\mathcal{X}`` and
``\\mathcal{X}``, we define the following discrete joint entropies:

- [`JointEntropyShannon`](@ref):
    ``H^S(X, Y) = -\\sum_{x\\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y) \\log p(x, y)``
    (Cover & Thomas, 2006)[^CoverThomas2006].
- [`JointEntropyRenyi`](@ref):
    ``H_q^R(X, Y) = -\\dfrac{1}{1-\\alpha} \\log \\sum_{i = 1}^N p_i^q``
    (Golshani et al., 2009)[^Golshani2009].
- [`JointEntropyTsallis`](@ref):
    ``H_q^T(X, Y) = -\\sum_{x\\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y)^q \\log_q p(x, y)``
    (Furuichi, 2006)[^Furuichi2006],
    where ``log_q(x, q) = \\dfrac{x^{1-q} - 1}{1-q}`` is the q-logarithm.

[^CoverThomas2006]:
    Thomas M. Cover and Joy A. Thomas. 2006. Elements of Information Theory (Wiley Series
    in Telecommunications and Signal Processing). Wiley-Interscience, USA.
[^Golshani2009]:
    Golshani, L., Pasha, E., & Yari, G. (2009). Some properties of Rényi
    entropy and Rényi entropy rate. Information Sciences, 179(14), 2426-2433.
[^Furuichi2006]:
    Furuichi, S. (2006). Information theoretical properties of Tsallis entropies. Journal
    of Mathematical Physics, 47(2), 023302.
"""
function entropy_joint(e::EntropyDefinition, args...)
    throw(ArgumentError("Joint entropy not defined and/or implemented for $e with $(args)"))
end

################################################################
# Types of joint entropy. Each type of joint entropy is its
# own type, so it can be used as a "module" in other measures.
################################################################
""" The supertype of all joint entropy measures. """
abstract type JointEntropy end

"""
    JointEntropyShannon <: JointEntropy
    JointEntropyShannon(; base = 2)

The Shannon joint entropy measure. See docstring of [`entropy_joint`](@ref) for definition.
"""
struct JointEntropyShannon{E<:Shannon} <: JointEntropy
    e::E
    function JointEntropyShannon(; base = 2)
        e = Shannon(; base)
        new{typeof(e)}(e)
    end
end

"""
    JointEntropyTsallis <: JointEntropy
    JointEntropyTsallis(; base = 2, q = 1.5)

The Tsallis joint entropy measure. See docstring of [`entropy_joint`](@ref) for definition.
"""
struct JointEntropyTsallis{E<:Tsallis} <: JointEntropy
    e::E
    function JointEntropyTsallis(; base = 2, q = 1.5)
        e = Tsallis(; base, q)
        new{typeof(e)}(e)
    end
end

"""
    JointEntropyRenyi <: JointEntropy
    JointEntropyRenyi(; base = 2, q = 1.5)

The Tsallis joint entropy measure. See docstring of [`entropy_joint`](@ref) for definition.
"""
struct JointEntropyRenyi{E<:Renyi} <: JointEntropy
    e::E
    function JointEntropyRenyi(; base = 2, q = 1.5)
        e = Renyi(; base, q)
        new{typeof(e)}(e)
    end
end

################################################################
# Discrete implementations
################################################################
function entropy_joint(measure::JointEntropyShannon, x...)
    # Define p(x...) log p(x...) = 0 if p(x....) = 0; (Cover & Thomas, 2006)
    # We circumvent this definition by directly counting *occurring pairs*.
    # Any non-occurring pair then gets probability zero automatically.
    X = StateSpaceSet(x...)
    return entropy(measure.e, CountOccurrences(), X)
end

function entropy_joint(measure::JointEntropyShannon, est::DifferentialEntropyEstimator, x...)
    X = StateSpaceSet(x...)
    return entropy(measure.e, est, X)
end

function entropy_joint(measure::JointEntropyRenyi, x...)
    # Direct analogue of Shannon version,
    #Golshani, L., Pasha, E., & Yari, G. (2009). Some properties of Rényi entropy and Rényi entropy rate. Information Sciences, 179(14), 2426-2433.
    X = StateSpaceSet(x...)
    return entropy(measure.e, CountOccurrences(), X)
end

function entropy_joint(measure::JointEntropyTsallis, x...)
    X = StateSpaceSet(x...)
    return entropy(measure.e, CountOccurrences(), X)
end

function entropy_joint(measure::JointEntropyShannon, c::ContingencyMatrix{T, 2}) where {T}
    base = measure.e.base
    h = 0.0
    for pij in c
        h += pij * log(pij)
    end
    return (-h) / log(base, ℯ)
end


# ``H_q^T(X, Y) = -\\sum_{x\\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y)^q \\log_q p(x, y)``
# (Furuichi, 2006)[^Furuichi2006],
# where ``log_q(x, q) = \\dfrac{x^{1-q} - 1}{1-q}`` is the q-logarithm.

function entropy_joint(measure::JointEntropyTsallis, c::ContingencyMatrix{T, 2}) where {T}
    base = measure.e.base
    q = measure.e.q
    h = 0.0
    for pij in c
        h += pij^q * logq(pij, q)
    end
    return (-h) / log(base, ℯ)
end


function entropy_joint(measure::JointEntropyRenyi, c::ContingencyMatrix{T, 2}) where {T}
    base = measure.e.base
    q = measure.e.q
    h = 0.0
    for pij in c
        h += pij^q * logq(pij, q)
    end
    return (-h) / log(base, ℯ)
end
