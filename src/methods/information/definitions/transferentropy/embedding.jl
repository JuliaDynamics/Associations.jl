export EmbeddingTE

"""
    EmbeddingTE(; dS = 1, dT = 1, dTf = 1, dC = 1, τS = -1, τT = -1, ηTf = 1, τC = -1)
    EmbeddingTE(opt::OptimiseTraditional, s, t, [c])

`EmbeddingTE` provide embedding parameters for transfer entropy analysis
using either [`TEShannon`](@ref), [`TERenyi`](@ref), or in general any subtype
of [`TransferEntropy`](@ref).

The second method finds parameters using the ["traditional"](https://juliadynamics.github.io/DynamicalSystems.jl/dev/embedding/traditional/)
optimised embedding techniques from DynamicalSystems.jl

## Convention for generalized delay reconstruction

We use the following convention.
Let ``s(i)`` be time series for the source variable, ``t(i)`` be the time series for the
target variable and ``c(i)`` the time series for the conditional variable. To compute
transfer entropy, we need the following marginals:

```math
\\begin{aligned}
T^{+} &= \\{t(i+\\eta^1), t(i+\\eta^2), \\ldots, (t(i+\\eta^{d_{T^{+}}}) \\} \\\\
T^{-} &= \\{ (t(i+\\tau^0_{T}), t(i+\\tau^1_{T}), t(i+\\tau^2_{T}), \\ldots, t(t + \\tau^{d_{T} - 1}_{T})) \\} \\\\
S^{-} &= \\{ (s(i+\\tau^0_{S}), s(i+\\tau^1_{S}), s(i+\\tau^2_{S}), \\ldots, s(t + \\tau^{d_{S} - 1}_{S})) \\} \\\\
C^{-} &= \\{ (c(i+\\tau^0_{C}), c(i+\\tau^1_{C}), c(i+\\tau^2_{C}), \\ldots, c(t + \\tau^{d_{C} - 1}_{C})) \\}
\\end{aligned}
```

Depending on the application, the delay reconstruction lags
``\\tau^k_{T} \\leq 0``, ``\\tau^k_{S} \\leq 0``, and ``\\tau^k_{C} \\leq 0``
may be equally spaced, or non-equally spaced. The same applied to the prediction lag(s),
but typically only a only a single predictions lag ``\\eta^k`` is used (so that
``d_{T^{+}} = 1``).

For transfer entropy, traditionally at least one ``\\tau^k_{T}``, one ``\\tau^k_{S}`` and
one ``\\tau^k_{C}`` equals zero. This way, the ``T^{-}``, ``S^{-}`` and ``C^{-}`` marginals
always contains present/past states, while the ``\\mathcal T`` marginal contain future
states relative to the other marginals. However, this is not a strict requirement,
and modern approaches that searches for optimal embeddings can return embeddings without
the intantaneous lag.

Combined, we get the generalized delay reconstruction
``\\mathbb{E} = (T^{+}_{(d_{T^{+}})}, T^{-}_{(d_{T})}, S^{-}_{(d_{S})}, C^{-}_{(d_{C})})``.
Transfer entropy is then computed as

```math
\\begin{aligned}
TE_{S \\rightarrow T | C} = \\int_{\\mathbb{E}} P(T^{+}, T^-, S^-, C^-)
\\log_{b}{\\left(\\frac{P(T^{+} | T^-, S^-, C^-)}{P(T^{+} | T^-, C^-)}\\right)},
\\end{aligned}
```

or, if conditionals are not relevant,

```math
\\begin{aligned}
TE_{S \\rightarrow T} = \\int_{\\mathbb{E}} P(T^{+}, T^-, S^-)
\\log_{b}{\\left(\\frac{P(T^{+} | T^-, S^-)}{P(T^{+} | T^-)}\\right)},
\\end{aligned}
```

Here,

- ``T^{+}`` denotes the ``d_{T^{+}}``-dimensional set of vectors furnishing the future
    states of ``T`` (almost always equal to 1 in practical applications),
- ``T^{-}`` denotes the ``d_{T}``-dimensional set of vectors furnishing the past and
    present states of ``T``,
- ``S^{-}`` denotes the ``d_{S}``-dimensional set of vectors furnishing the past and
    present of ``S``, and
- ``C^{-}`` denotes the ``d_{C}``-dimensional set of vectors furnishing the past and
    present of ``C``.

## Keyword arguments

- `dS`, `dT`, `dC`, `dTf` (`f` for *future*) are the dimensions of the ``S^{-}``,
    ``T^{-}``, ``C^{-}`` and ``T^{+}`` marginals. The parameters `dS`, `dT`, `dC` and `dTf`
    must each be a *positive* integer number.
-  `τS`, `τT`, `τC` are the embedding lags for ``S^{-}``, ``T^{-}``, ``C^{-}``.
    Each parameter are integers `∈ 𝒩⁰⁻`, or a vector of integers `∈ 𝒩⁰⁻`, so
    that ``S^{-}``, ``T^{-}``, ``C^{-}`` always represents present/past values.
    If e.g. `τT` is an integer, then for the ``T^-`` marginal is constructed using
    lags ``\\tau_{T} = \\{0, \\tau, 2\\tau, \\ldots, (d_{T}- 1)\\tau_T \\}``.
    If is a vector, e.g. `τΤ = [-1, -5, -7]`, then the dimension `dT` must match the lags,
    and precisely those lags are used: ``\\tau_{T} = \\{-1, -5, -7 \\}``.
- The prediction lag(s) `ηTf` is a positive integer. Combined with the requirement
    that the other delay parameters are zero or negative, this ensures that we're
    always predicting from past/present to future. In typical applications,
    `ηTf = 1` is used for transfer entropy.

## Examples

Say we wanted to compute the Shannon transfer entropy
``TE^S(S \\to T) = I^S(T^+; S^- | T^-)``. Using some modern procedure for
determining optimal embedding parameters using
[methods from DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/dev/embedding/unified/),
we find that the optimal embedding of ``T^{-}`` is three-dimensional and is given by
the lags `[0, -5, -8]`. Using the same procedure, we find that the optimal embedding
of ``S^{-}`` is two-dimensional with lags ``[-1, -8]``. We want to predicting a univariate
version of the target variable one time step into the future (`ηTf = 1`).
The total embedding is then the set of embedding vectors

``E_{TE} = \\{ (T(i+1), S(i-1), S(i-8), T(i), T(i-5), T(i-8)) \\}``. Translating
this to code, we get:

```jldoctest
using CausalityTools
julia> EmbeddingTE(dT=3, τT=[0, -5, -8], dS=2, τS=[-1, -4], ηTf=1)

# output
EmbeddingTE(dS=2, dT=3, dC=1, dTf=1, τS=[-1, -4], τT=[0, -5, -8], τC=-1, ηTf=1)
```
"""
@Base.kwdef struct EmbeddingTE
    dS::Union{Int, AbstractVector{Int}} = 1
    dT::Union{Int, AbstractVector{Int}} = 1
    dTf::Union{Int, AbstractVector{Int}} = 1
    dC::Union{Int, AbstractVector{Int}, Nothing} = 1
    τS::Union{Int, AbstractVector{Int}} = -1
    τT::Union{Int, AbstractVector{Int}} = -1
    ηTf::Union{Int, AbstractVector{Int}} = 1
    τC::Union{Int, AbstractVector{Int}, Nothing} = -1

    function EmbeddingTE(
            dS::Union{Int, AbstractVector{Int}},
            dT::Union{Int, AbstractVector{Int}},
            dTf::Union{Int, AbstractVector{Int}},
            dC::Union{Int, AbstractVector{Int}},
            τS::Union{Int, AbstractVector{Int}},
            τT::Union{Int, AbstractVector{Int}},
            ηTf::Union{Int, AbstractVector{Int}},
            τC::Union{Int, AbstractVector{Int}})

        if dS isa Int
            dS > 0 || throw(ArgumentError("dimension for marginal S must be a positive integer (got dS=$(dS))"))
        end
        if dT isa Int
            dT > 0 || throw(ArgumentError("dimension for marginal T must be a positive integer (got dT=$(dT))"))
        end
        if dC isa Int
            dC > 0 || throw(ArgumentError("dimension for marginal C must be a positive integer (got dC=$(dC))"))
        end
        if dTf isa Int
            dTf > 0 || throw(ArgumentError("dimension for marginal Tf must be a positive integer (got dTf=$(dTf))"))
        end
        if τS isa Int
            τS < 0 || throw(ArgumentError("delay for marginal S must be a negative integer (got τS=$(τS))"))
        end
        if τT isa Int
            τT < 0 || throw(ArgumentError("delay for marginal T must be a negative integer (got τT=$(τT))"))
        end
        if τC isa Int
            τC < 0 || throw(ArgumentError("delay for marginal C must be a negative integer (got τC=$(τC))"))
        end

        if τS isa AbstractVector{Int} || τS isa AbstractUnitRange{Int64}
            all(τS .<= 0) || throw(ArgumentError("delays for marginal S must be <= 0 (got τS=$(τS))"))
        end

        if τT isa AbstractVector{Int} || τT isa AbstractUnitRange{Int64}
            all(τT .<= 0) || throw(ArgumentError("delays for marginal T must be <= 0 (got τT=$(τT))"))
        end

        if τC isa AbstractVector{Int} || τC isa AbstractUnitRange{Int64}
            all(τC .<= 0) || throw(ArgumentError("delays for marginal C must be <= 0 (got τC=$(τC))"))
        end

        new(dS, dT, dTf, dC, τS, τT, ηTf, τC)
    end

end

function Base.show(io::IO, x::EmbeddingTE)
    s = "EmbeddingTE(dS=$(x.dS), dT=$(x.dT), dC=$(x.dC), dTf=$(x.dTf), τS=$(x.τS), τT=$(x.τT), τC=$(x.τC), ηTf=$(x.ηTf))"
    print(io, s)
end

#include("optimization/traditional_optimal_embedding.jl")