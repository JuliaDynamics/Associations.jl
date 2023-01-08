export EmbeddingTE

"""
    EmbeddingTE(; dS = 1, dT = 1, d𝒯 = 1, dC = 1, τS = -1, τT = -1, η𝒯 = 1, τC = -1)

`EmbeddingTE` provide embedding parameters for transfer entropy analysis
using either [`TEShannon`](@ref), [`TERenyi`](@ref), or in general any subtype
of [`TransferEntropy`](@ref), which in turns dictates the embedding used with
[`transferentropy`](@ref).

## Convention for generalized delay reconstruction

We use the following convention.
Let ``x(t)`` be time series for the source variable, ``y(t)`` be the time series for the
target variable and ``z(t)`` the time series for any conditional variable. To compute
transfer entropy, we need the following marginals:

```math
\\begin{aligned}
\\mathcal{T}^{(d_{\\mathcal{T}})} &= \\{(y(t+\\eta^{d_{\\mathcal{T}}}), \\ldots, y(t+\\eta^2), y(t+\\eta^1) \\} \\\\
T^{(d_{T})} &= \\{ (y(t+\\tau^0_{T}), y(t+\\tau^1_{T}), y(t+\\tau^2_{T}), \\ldots, y(t + \\tau^{d_{T} - 1}_{T})) \\} \\\\
S^{(d_{S})} &= \\{ (x(t+\\tau^0_{S}), x(t+\\tau^1_{S}), x(t+\\tau^2_{S}), \\ldots, x(t + \\tau^{d_{S} - 1}_{S})) \\} \\\\
C^{(d_{C})} &= \\{ (z(t+\\tau^0_{C}), z(t+\\tau^1_{C}), z(t+\\tau^2_{C}), \\ldots, z(t + \\tau^{d_{C} - 1}_{C})) \\}
\\end{aligned}
```

Depending on the application, the delay reconstruction lags ``\\tau^k_{T} \\leq 0``, ``\\tau^k_{S} \\leq 0``, and ``\\tau^k_{C} \\leq 0``
may be equally spaced, or non-equally spaced. The predictions lags ``\\eta^k``may also be equally spaced
or non-equally spaced, but are always positive. For transfer entropy, convention dictates that at least one
``\\tau^k_{T}``, one ``\\tau^k_{S}`` and one ``\\tau^k_{C}`` equals zero. This way, the ``T``, ``S`` and ``C`` marginals
always contains present/past states,
while the ``\\mathcal T`` marginal contain future states relative to the other marginals.

Combined, we get the generalized delay reconstruction ``\\mathbb{E} = (\\mathcal{T}^{(d_{\\mathcal{T}})}, T^{(d_{T})}, S^{(d_{S})}, C^{(d_{C})})``. Transfer entropy is then computed as

```math
\\begin{aligned}
TE_{S \\rightarrow T | C} = \\int_{\\mathbb{E}} P(\\mathcal{T}, T, S, C) \\log_{b}{\\left(\\frac{P(\\mathcal{T} | T, S, C)}{P(\\mathcal{T} | T, C)}\\right)},
\\end{aligned}
```

or, if conditionals are not relevant,

```math
\\begin{aligned}
TE_{S \\rightarrow T} = \\int_{\\mathbb{E}} P(\\mathcal{T}, T, S) \\log_{b}{\\left(\\frac{P(\\mathcal{T} | T, S)}{P(\\mathcal{T} | T)}\\right)},
\\end{aligned}
```

Here,

- ``\\mathcal{T}`` denotes the ``d_{\\mathcal{T}}``-dimensional set of vectors furnishing the future states of ``T``,
- ``T`` denotes the ``d_{T}``-dimensional set of vectors furnishing the past and present states of ``T``,
- ``S`` denotes the ``d_{S}``-dimensional set of vectors furnishing the past and present of ``S``, and
- ``C`` denotes the ``d_{C}``-dimensional set of vectors furnishing the past and present of ``C``.

## Keyword arguments

### Specifying dimensions for generalized delay reconstructions of marginals

`dS`, `dT`, `d𝒯`, and `dC` are the dimensions of the ``S``, ``T``, ``\\mathcal{T}``,
and ``C`` marginals. The dimensions of each marginal can be specified manually by setting
either `dS`, `dT`, `d𝒯`, or `dC` to a *positive* integer number. Alternatively, the dimension
of each marginal can be optimised by setting either `dS`, `dT`, `d𝒯`, or `dC` to an
instance of [`OptimiseDim`](@ref)
(e.g. `EmbeddingTE(dT = OptimDim(method_delay = "ac_zero", method_dim = "f1nn")`).

### Specifying delays for generalized delay reconstructions of marginals

The corresponding embedding delay lags are given by `τS`, `τT` and `τC`. The delays
for each marginal can be specified manually by setting either `dS`, `dT`, `d𝒯`, or `dC`
to a *negative* integer number. The delay defaults for each marginal is -1 (but is set to zero
if the marginal is one-dimensional), and must always be negative. Alternatively, delays can
be estimated numerically by setting either `dS`, `dT`, `d𝒯`, and `dC`
to an instance of [`OptimiseDelay`](@ref) (e.g. `dS = OptimiseDelay(method_delay = "ac_zero")`).

The prediction lag `η` can be either positive or negative, but should not be zero.

In summary, one can provide

- A single delay ``\\tau``, in which case ``\\tau_{T} = \\{0, \\tau, 2\\tau, \\ldots, (d_{T}- 1)\\tau \\}``, or
- All the delays manually. If so, then the number of delays must match the dimension of the marginal).

For the prediction lag, one can provide

- A single delay ``\\eta_f``, in which case ``\\eta_{\\mathcal{T}} = \\{\\eta_f, 2\\eta_f, \\ldots, (d_{\\mathcal{T}} - 1)\\eta_f \\}``, or
- All the delays manually. If so, then the number of delays must equal ``d_{\\mathcal{T}}``, which is the dimension of the marginal).

!!! note
    If both the delay and the dimension for a given marginal is to be estimated numerically, make sure
    to use the same delay estimation method for both
    the [`OptimiseDelay`](@ref) and  [`OptimiseDim`](@ref) instances.

## Examples

Say we wanted to compute the Shannon transfer entropy
``TE^S{S \\to T} = I^S(T^+; S^- | T^-)``. Using some modern procedure for
determining optimal embedding parameters using
[methods from DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/dev/embedding/unified/),
we find that the optimal embedding of ``T`` is three-dimensional and is given by
the lags `[0, -5, -8]`. Using the same procedure, we find that the optimal embedding
of ``S`` is two-dimensional with lags ``[-1, -8]``. We want to predicting a univariate
version of the target variable one time step into the future (`η𝒯 = 1`).
The total embedding is then the set of embedding vectors

``E_{TE} = \\{ (T(i+1), S(i-1), S(i-8), T(i), T(i-5), T(i-8)) \\}``. Translating
this to code, we get:

```jldoctest
using CausalityTools
julia> EmbeddingTE(dT=3, τT=[0, -5, -8], dS=2, τS=[-1, -4], η𝒯=1)

# output
EmbeddingTE(dS=2, dT=3, dC=1, d𝒯=1, τS=[-1, -4], τT=[0, -5, -8], τC=-1, η𝒯=1)
```
"""
@Base.kwdef struct EmbeddingTE
    dS::Union{Int, AbstractVector{Int}} = 1
    dT::Union{Int, AbstractVector{Int}} = 1
    d𝒯::Union{Int, AbstractVector{Int}} = 1
    dC::Union{Int, AbstractVector{Int}, Nothing} = 1
    τS::Union{Int, AbstractVector{Int}} = -1
    τT::Union{Int, AbstractVector{Int}} = -1
    η𝒯::Union{Int, AbstractVector{Int}} = 1
    τC::Union{Int, AbstractVector{Int}, Nothing} = -1

    function EmbeddingTE(
            dS::Union{Int, AbstractVector{Int}},
            dT::Union{Int, AbstractVector{Int}},
            d𝒯::Union{Int, AbstractVector{Int}},
            dC::Union{Int, AbstractVector{Int}},
            τS::Union{Int, AbstractVector{Int}},
            τT::Union{Int, AbstractVector{Int}},
            η𝒯::Union{Int, AbstractVector{Int}},
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
        if d𝒯 isa Int
            d𝒯 > 0 || throw(ArgumentError("dimension for marginal 𝒯 must be a positive integer (got d𝒯=$(d𝒯))"))
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

        new(dS, dT, d𝒯, dC, τS, τT, η𝒯, τC)
    end

end

function Base.show(io::IO, x::EmbeddingTE)
    s = "EmbeddingTE(dS=$(x.dS), dT=$(x.dT), dC=$(x.dC), d𝒯=$(x.d𝒯), τS=$(x.τS), τT=$(x.τT), τC=$(x.τC), η𝒯=$(x.η𝒯))"
    print(io, s)
end
