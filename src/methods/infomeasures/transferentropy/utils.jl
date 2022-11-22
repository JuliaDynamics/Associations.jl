import DelayEmbeddings: genembed, AbstractDataset, Dataset
export EmbeddingTE

function rc(x::Union{AbstractDataset, AbstractVector{T}},
        dim::Union{Int, AbstractVector{Int}},
        τ::Union{Int, AbstractVector{Int}}, forward = false) where T <: Union{Number, AbstractVector}

    if typeof(x) <: AbstractVector{T} where T <: Number
        # Fixed dimension, fixed lag
        if dim isa Int && τ isa Int
            pos = [1 for x in dim:-1:1]
            if τ > 0
                lags = [d*τ for d in dim:-1:1]
            elseif τ <= 0 && forward
                lags = [d*τ for d in dim:-1:1]
            elseif τ <= 0 && !forward
                lags = [d*τ for d in 0:(dim-1)]
            end
        end

        # Fixed dimension, multiple lags (number of lags must match dimension)
        if dim isa Int && length(τ) > 1
            length(τ) == dim || throw(ArgumentError("length(τ) must equal dim if multiple lags are specified manually (got length(τ) = $(length(τ)), dim=$(dim))"))
            length(unique(τ)) == length(τ) || throw(ArgumentError("Tried to specify reconstruction lags manually, but there are repeated lags. Provide unique lags."))
            pos = pos = [1 for x in dim:-1:1]
            return pos, τ
        end
    end

    # Multiple time series input
    if (x isa AbstractVector{T} where T <: AbstractVector{N} where N <: Number) || (x isa AbstractDataset)
        if x isa AbstractDataset
            N = length(x)
        elseif x isa AbstractVector
            N = size(x, 1)
        end


        if dim isa Int
            dim % N == 0 || throw(ArgumentError("If using multiple (`N` different) time series in a marginal, each time series is lagged `dim/N` times. Hence, `dim` must be a multiple of `N`."))

            lD = Int(dim / N) # "local" reconstruction dimension for this time series

            pos = Vector{Int}(undef, 0)

            for i = 1:N
                append!(pos, repeat([i], lD))
            end

            if τ isa Int
                if τ > 0
                    llags = [d*τ for d in lD:-1:1]
                elseif τ <= 0 && forward
                    llags = [d*τ for d in lD:-1:1]
                elseif τ <= 0 && !forward
                    llags = [d*τ for d in 0:(lD - 1)]
                end
                lags = repeat(llags, Int(dim/lD))
            elseif τ isa AbstractVector{Int}
                length(τ) == N || throw(ArgumentError("Tried using `N = $N` different time series in a marginal, but $(length(τ)) reconstruction lags were provided. The number of lags must equal `N`"))
                lags = Vector{Int}(undef, 0)
                for i = 1:length(x)
                    if τ[i] > 0
                        llags = [d*τ[i] for d in lD:-1:1]
                    elseif τ[i] <= 0 && forward
                        llags = [d*τ[i] for d in 0:(lD - 1)]
                    elseif τ[i] <= 0 && !forward
                        llags = [d*τ[i] for d in 0:(lD - 1)]
                    end
                    append!(lags, llags)
                end
            end
        end

        if dim isa AbstractVector{Int}
            length(dim) == N || throw(ArgumentError("There must be precisely one dimension specification per time series. Got $(length(dim)) specifications for $N time series."))
            pos = Vector{Int}(undef, 0)
            lags = Vector{Int}(undef, 0)

            for (i, lD) in enumerate(dim)
                append!(pos, repeat([i], lD))

                if τ isa Int
                    if τ > 0
                        llags = [d*τ for d in lD:-1:1]
                    elseif τ <= 0 && !forward
                        llags = [d*τ for d in lD:-1:1]
                    elseif τ <= 0 && forward
                        llags = [d*τ for d in 0:(lD - 1)]
                    end
                    append!(lags, llags)
                elseif τ isa AbstractVector{Int}
                    length(τ) == N || throw(ArgumentError("Tried using `N = $N` different time series in a marginal, but $(length(τ)) reconstruction lags were provided. The number of lags must equal `N`"))
                    if τ[i] > 0
                        llags = [d*τ[i] for d in lD:-1:1]
                    elseif τ[i] <= 0 && forward
                        llags = [d*τ[i] for d in lD:-1:1]
                    elseif τ[i] <= 0 && !forward
                        llags = [d*τ[i] for d in 0:(lD - 1)]
                    end
                    append!(lags, llags)
                end
            end
        end
    end
    return pos, lags
end


"""
    EmbeddingTE(; dS = 1, dT = 1, d𝒯 = 1, dC = 1, τS = -1, τT = -1, η𝒯 = 1, τC = -1)

Embedding parameters for transfer entropy analysis.

## Convention for generalized delay reconstruction

This struct contains instructions for transfer entropy computations using the following convention.
Let ``x(t)`` be time series for the source variable, ``y(t)`` be the time series for the target variable and
``z(t)`` the time series for any conditional variable. To compute transfer entropy, we need the
following marginals:


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

Default parameters:

```jldoctest
using TransferEntropy
p = EmbeddingTE()

# output
EmbeddingTE(dS=1, dT=1, dC=1, d𝒯=1, τS=-1, τT=-1, τC=-1, η𝒯=1)
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

function get_delay_reconstruction_params(source, target, p::EmbeddingTE)
    pos_𝒯, lags_𝒯 = rc(target, p.d𝒯, p.η𝒯, true)
    pos_T, lags_T = rc(target, p.dT, p.τT, false)
    pos_S, lags_S = rc(source, p.dS, p.τS, false)
    pos_C, lags_C = rc(source, p.dC, p.τC, false)

    js = ([pos_𝒯; pos_T; pos_S]...,)
    τs = ([lags_𝒯; lags_T; lags_S]...,)

    return τs, js
end

function get_delay_reconstruction_params(source, target, cond, p::EmbeddingTE)
    pos_𝒯, lags_𝒯 = rc(target, p.d𝒯, p.η𝒯, true)
    pos_T, lags_T = rc(target, p.dT, p.τT, false)
    pos_S, lags_S = rc(source, p.dS, p.τS, false)
    pos_C, lags_C = rc(cond, p.dC, p.τC, false)

    js = ([pos_𝒯; pos_T; pos_S; pos_C]...,)
    τs = ([lags_𝒯; lags_T; lags_S; pos_C]...,)

    return τs, js
end

"""
    te_embed(source::AbstractVector{T}, target::AbstractVector{T}, p::EmbeddingTE) → (points, vars, τs)
    te_embed(source::AbstractVector{T}, target::AbstractVector{T}, cond::AbstractVector{T}, p::EmbeddingTE) → (points, vars, τs)

Generalised delay reconstruction of `source` and `target` (and `cond` if provided)
for transfer entropy computation using embedding parameters provided by the [`EmbeddingTE`](@ref)
instance `p`.

Returns a tuple of the embedded `points`, `vars` (a [`TEVars`](@ref) instance that keeps track of which
variables of the embedding belong to which marginals of the reconstruction; indices are: source = 1,
target = 2, cond = 3), and a tuple `τs`, which stores the lags for each variable of the reconstruction.
"""
function te_embed(source::AbstractVector{T}, target::AbstractVector{T}, p::EmbeddingTE) where T

    #@show p.τS
    #if (p.τS isa Int && p.τS > 0) || (length(p.τS) > 1 && any(p.τS[p.τS .> 0]))
    #    @warn("Backwards lag τS should be negative. You might be getting nonsensical results!")
    #end

    # Get lags and posisions separately for each marginal
    pos_𝒯, lags_𝒯 = rc(target, p.d𝒯, p.η𝒯, true)
    pos_T, lags_T = rc(target, p.dT, p.τT, false)
    pos_S, lags_S = rc(source, p.dS, p.τS, false)

    # Add one to the index of the positions for the target (rc doesn't know it is in fact our second time series)
    # TODO: make sure this works when `source` and `target` are multiple time series
    pos_𝒯 .= pos_𝒯 .+ 1
    pos_T .= pos_T .+ 1

    js = ([pos_𝒯; pos_T; pos_S]...,)
    τs = ([lags_𝒯; lags_T; lags_S]...,)

    # TODO: This only works for single time series at the moment
    ts = Dataset(source, target)

    # The reconstructed points
    pts = genembed(ts, τs, js)
    d𝒯 = length(pos_𝒯)
    dT = length(pos_T)
    dS = length(pos_S)

    # Which columns/variables map to which marginals?
    vars = TEVars(
        𝒯  = 1:(d𝒯)           |> collect,
        T = 1+(d𝒯):dT+(d𝒯)     |> collect,
        S = 1+(dT+d𝒯):dS+(d𝒯+dT) |> collect)

    return pts, vars, τs, js
end

function te_embed(source::AbstractVector{T}, target::AbstractVector{T}, cond::AbstractVector{T}, p::EmbeddingTE) where T

    #@show p.τS
    #if (p.τS isa Int && p.τS > 0) || (length(p.τS) > 1 && any(p.τS[p.τS .> 0]))
    #    @warn("Backwards lag τS should be negative. You might be getting nonsensical results!")
    #end
    # Get lags and posisions separately for each marginal
    pos_𝒯, lags_𝒯 = rc(target, p.d𝒯, p.η𝒯,  true)
    pos_T, lags_T = rc(target, p.dT, p.τT, false)
    pos_S, lags_S = rc(source, p.dS, p.τS, false)
    pos_C, lags_C = rc(cond,   p.dC, p.τC, false)

    # Add one to the index of the positions for the target (rc doesn't know it is in fact our second time series)
    # TODO: make sure this works when `source` and `target` are multiple time series
    pos_𝒯 .= pos_𝒯 .+ 1
    pos_T .= pos_T .+ 1
    pos_C .= pos_C .+ 2

    js = ([pos_𝒯; pos_T; pos_S; pos_C]...,)
    τs = ([lags_𝒯; lags_T; lags_S; lags_C]...,)

    # TODO: This only works for single time series at the moment
    ts = Dataset(source, target, cond)

    # The reconstructed points
    pts = genembed(ts, τs, js)
    d𝒯 = length(pos_𝒯)
    dT = length(pos_T)
    dS = length(pos_S)
    dC = length(pos_C)

    # Which columns/variables map to which marginals?
    vars = TEVars(
        𝒯 = 1:(d𝒯)               |> collect,
        T = 1+(d𝒯):dT+(d𝒯)         |> collect,
        S = 1+(dT+d𝒯):dS+(d𝒯+dT)     |> collect,
        C = 1+(dT+d𝒯+dS):dC+(d𝒯+dT+dS) |> collect)

    return pts, vars, τs, js
end


"""
    TEVars(𝒯::Vector{Int}, T::Vector{Int}, S::Vector{Int})
    TEVars(𝒯::Vector{Int}, T::Vector{Int}, S::Vector{Int}, C::Vector{Int})
    TEVars(;𝒯 = Int[], T = Int[], S = Int[], C = Int[]) -> TEVars

Which axes of the state space correspond to the future of the target (`𝒯`),
the present/past of the target (`T`), the present/past of the source (`S`), and
the present/past of any conditioned variables (`C`)?  This information is used by
the transfer entropy estimators to ensure that marginal distributions are computed correctly.

Indices correspond to variables of the embedding, or, equivalently, colums of a `Dataset`.

- The three-argument constructor assumes there will be no conditional variables.
- The four-argument constructor assumes there will be conditional variables.

"""
struct TEVars
    𝒯::Vector{Int}
    T::Vector{Int}
    S::Vector{Int}
    C::Vector{Int}
end

"""
    TEVars(𝒯::Vector{Int},
            T::Vector{Int},
            S::Vector{Int})

Which axes of the state space correspond to the future of the target,
the present/past of the target, and the present/past of the source?
Indices correspond to variables of the embedding or colums of a `Dataset`.

This information is used by the transfer entropy estimators to ensure
the marginal distributions are computed correctly.

This three-argument constructor assumes there will be no conditional variables.
"""
TEVars(𝒯::Vector{Int}, T::Vector{Int}, S::Vector{Int}) = TEVars(𝒯, T, S, Int[])

TEVars(;𝒯::Vector{Int} = Int[],
	    T::Vector{Int} = Int[],
	    S::Vector{Int} = Int[],
		C::Vector{Int} = Int[]) =
	TEVars(𝒯, T, S, C)

function Base.show(io::IO, tv::TEVars)
    s = "$(typeof(tv))(𝒯 = $(tv.𝒯), T = $(tv.T), S = $(tv.S), C = $(tv.C))"
    print(io, s)
end

function get_marginals(measure::TransferEntropy, s, t; emb::EmbeddingTE)
    pts, vars, τs, js = te_embed(s, t, emb)

    # Get marginals
    ST = pts[:, [vars.S; vars.T]]
    T𝒯 = pts[:, [vars.𝒯; vars.T]]
    T = pts[:, vars.T]
    joint = pts

    return joint, ST, T𝒯, T
end

function get_marginals(measure::TransferEntropy, s, t, c; emb::EmbeddingTE)
    pts, vars, τs, js = te_embed(s, t, c, emb)

    # Get marginals
    ST = pts[:, [vars.S; vars.T; vars.C]]
    T𝒯 = pts[:, [vars.𝒯; vars.T; vars.C]]
    T = pts[:, [vars.T; vars.C]]
    joint = pts

    return joint, ST, T𝒯, T
end

# map a set of pre-embedded points to the correct marginals for transfer entropy computation
function get_marginals(measure::TransferEntropy, pts::AbstractDataset; emb::TEVars)
    # Get marginals
    ST = pts[:, [vars.S; vars.T]]
    T𝒯 = pts[:, [vars.𝒯; vars.T]]
    T = pts[:, vars.T]
    joint = pts

    return joint, ST, T𝒯, T
end
