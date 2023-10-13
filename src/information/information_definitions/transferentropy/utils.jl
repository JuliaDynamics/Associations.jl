import DelayEmbeddings: genembed, AbstractStateSpaceSet, StateSpaceSet
export EmbeddingTE

function rc(x::Union{AbstractStateSpaceSet, AbstractVector{T}},
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
    if (x isa AbstractVector{T} where T <: AbstractVector{N} where N <: Number) || (x isa AbstractStateSpaceSet)
        if x isa AbstractStateSpaceSet
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


function get_delay_reconstruction_params(source, target, p::EmbeddingTE)
    pos_Tf, lags_Tf = rc(target, p.dTf, p.ηTf, true)
    pos_T, lags_T = rc(target, p.dT, p.τT, false)
    pos_S, lags_S = rc(source, p.dS, p.τS, false)
    pos_C, lags_C = rc(source, p.dC, p.τC, false)

    js = ([pos_Tf; pos_T; pos_S]...,)
    τs = ([lags_Tf; lags_T; lags_S]...,)

    return τs, js
end

function get_delay_reconstruction_params(source, target, cond, p::EmbeddingTE)
    pos_Tf, lags_Tf = rc(target, p.dTf, p.ηTf, true)
    pos_T, lags_T = rc(target, p.dT, p.τT, false)
    pos_S, lags_S = rc(source, p.dS, p.τS, false)
    pos_C, lags_C = rc(cond, p.dC, p.τC, false)

    js = ([pos_Tf; pos_T; pos_S; pos_C]...,)
    τs = ([lags_Tf; lags_T; lags_S; pos_C]...,)

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
function te_embed(p::EmbeddingTE, source::AbstractVector{T}, target::AbstractVector{T}) where T

    #@show p.τS
    #if (p.τS isa Int && p.τS > 0) || (length(p.τS) > 1 && any(p.τS[p.τS .> 0]))
    #    @warn("Backwards lag τS should be negative. You might be getting nonsensical results!")
    #end

    # Get lags and posisions separately for each marginal
    pos_Tf, lags_Tf = rc(target, p.dTf, p.ηTf, true)
    pos_T, lags_T = rc(target, p.dT, p.τT, false)
    pos_S, lags_S = rc(source, p.dS, p.τS, false)

    # Add one to the index of the positions for the target (rc doesn't know it is in fact our second time series)
    # TODO: make sure this works when `source` and `target` are multiple time series
    pos_Tf .= pos_Tf .+ 1
    pos_T .= pos_T .+ 1

    js = ([pos_Tf; pos_T; pos_S]...,)
    τs = ([lags_Tf; lags_T; lags_S]...,)

    # TODO: This only works for single time series at the moment
    ts = StateSpaceSet(source, target)

    # The reconstructed points
    pts = genembed(ts, τs, js)
    dTf = length(pos_Tf)
    dT = length(pos_T)
    dS = length(pos_S)

    # Which columns/variables map to which marginals?
    vars = TEVars(
        Tf  = 1:(dTf)           |> collect,
        T = 1+(dTf):dT+(dTf)     |> collect,
        S = 1+(dT+dTf):dS+(dTf+dT) |> collect)

    return pts, vars, τs, js
end

function te_embed(p::EmbeddingTE, source::AbstractVector{T}, target::AbstractVector{T}, cond::AbstractVector{T}) where T

    #@show p.τS
    #if (p.τS isa Int && p.τS > 0) || (length(p.τS) > 1 && any(p.τS[p.τS .> 0]))
    #    @warn("Backwards lag τS should be negative. You might be getting nonsensical results!")
    #end
    # Get lags and posisions separately for each marginal
    pos_Tf, lags_Tf = rc(target, p.dTf, p.ηTf,  true)
    pos_T, lags_T = rc(target, p.dT, p.τT, false)
    pos_S, lags_S = rc(source, p.dS, p.τS, false)
    pos_C, lags_C = rc(cond,   p.dC, p.τC, false)

    # Add one to the index of the positions for the target (rc doesn't know it is in fact our second time series)
    # TODO: make sure this works when `source` and `target` are multiple time series
    pos_Tf .= pos_Tf .+ 1
    pos_T .= pos_T .+ 1
    pos_C .= pos_C .+ 2

    js = ([pos_Tf; pos_T; pos_S; pos_C]...,)
    τs = ([lags_Tf; lags_T; lags_S; lags_C]...,)

    # TODO: This only works for single time series at the moment
    ts = StateSpaceSet(source, target, cond)

    # The reconstructed points
    pts = genembed(ts, τs, js)
    dTf = length(pos_Tf)
    dT = length(pos_T)
    dS = length(pos_S)
    dC = length(pos_C)

    # Which columns/variables map to which marginals?
    vars = TEVars(
        Tf = 1:(dTf)               |> collect,
        T = 1+(dTf):dT+(dTf)         |> collect,
        S = 1+(dT+dTf):dS+(dTf+dT)     |> collect,
        C = 1+(dT+dTf+dS):dC+(dTf+dT+dS) |> collect)

    return pts, vars, τs, js
end

"""
    TEVars(Tf::Vector{Int}, T::Vector{Int}, S::Vector{Int})
    TEVars(Tf::Vector{Int}, T::Vector{Int}, S::Vector{Int}, C::Vector{Int})
    TEVars(;Tf = Int[], T = Int[], S = Int[], C = Int[]) -> TEVars

Which axes of the state space correspond to the future of the target (`Tf`),
the present/past of the target (`T`), the present/past of the source (`S`), and
the present/past of any conditioned variables (`C`)?  This information is used by
the transfer entropy estimators to ensure that marginal distributions are computed correctly.

Indices correspond to variables of the embedding, or, equivalently, colums of a `StateSpaceSet`.

- The three-argument constructor assumes there will be no conditional variables.
- The four-argument constructor assumes there will be conditional variables.

"""
Base.@kwdef struct TEVars
    Tf::Vector{Int} = Int[]
    T::Vector{Int} = Int[]
    S::Vector{Int} = Int[]
    C::Vector{Int} = Int[]
    TEVars(𝒯, T, S, C) = new(𝒯, T, S, C)
    TEVars(𝒯, T, S) = new(𝒯, T, S, Int[])
end

function Base.show(io::IO, tv::TEVars)
    s = "$(typeof(tv))(Tf = $(tv.Tf), T = $(tv.T), S = $(tv.S), C = $(tv.C))"
    print(io, s)
end
