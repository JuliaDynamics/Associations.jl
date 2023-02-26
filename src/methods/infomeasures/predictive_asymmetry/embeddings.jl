using DelayEmbeddings
export Pecuzal
export optimize_embedding

"""
    Pecuzal(; dmax = 3, τmax = 30, w = 0)

Used with [`optimize_embedding`](@ref) to find the optimal embedding lags using
the Pecuzal algorithm. `w` is the Theiler window. `dmax` is the maximum
dimension and we query lags `1:τmax`.
"""
Base.@kwdef struct Pecuzal
    dmax::Int = 3
    τmax::Int = 30
    w::Int = 0
end


function lags_pecuzal(x::Union{Vector, AbstractDataset{1}}, Tmax::Int = 30;
        w = 0)

    Y_d, τ_vals_d, ts_vals_d, Ls_d , ε★_d = pecuzal_embedding(x; τs = 0:Tmax, w,
        verbose = false)
    # If valid embedding is not achieved, default to embedding lag 1.
    if length(τ_vals_d) == 1 && τ_vals_d[1] == 0
        return 1
    else
        # We only want strictly positive lags.
        return τ_vals_d[τ_vals_d .> 0]
    end
end

function optimize_embedding(method::Pecuzal, s::AbstractVector)
    rS = estimate_delay(s, "mi_min", 1:method.τmax)
    τS = unique(lags_pecuzal(s, max(rS, method.τmax)))
    τS = first(τS[τS .> 0], method.dmax)
    return τS
end
function optimize_embedding(method::Pecuzal, s::AbstractDataset{1})
    τS = unique(lags_pecuzal(s, method.τmax))
    τS = first(τS[τS .> 0], method.dmax)
    return τS
end

export optimize_embeddings
