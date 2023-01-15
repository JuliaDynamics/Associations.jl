using DelayEmbeddings

export OptimiseTraditional
export optimize_marginals

"""
    OptimiseTraditional(; η = 1, maxlag = 50, maxdim = 10,
        method = delay_ifnn, dmethod = "mi_min")

Optimize embedding parameters using traditional delay embedding optimization methods
with a maximum lag of `maxlag` and a maximum dimension of `maxdim`. `method` can
be either `delay_ifnn`, `delay_fnn` or `delay_f1nn`.
"""
Base.@kwdef struct OptimiseTraditional{L}
    η::Int = 1
    maxlag::L = 0.05
    maxdim::Int = 5
    method::Function = delay_ifnn
    dmethod::AbstractString = "mi_min"
end

"""
    optimize_marginals([scheme = OptimiseTraditional()], s, t, [c])

Optimize marginal embeddings for transfer entropy computation from source time series `s`
to target time series `t`, conditioned on `c` if `c` is given, using the provided
optimization `scheme`.
"""
function optimize_marginals end
getlags(opt::OptimiseTraditional{<:Float64}, x) = 1:floor(Int, length(x)*opt.maxlag)
getlags(opt::OptimiseTraditional{<:Int}, args...) = 1:opt.maxlag

function optimize_marginals(opt::OptimiseTraditional, s, t)
    τs = getlags(opt, t)
    dims = 1:opt.maxdim
    τT = estimate_delay(t, opt.dmethod, τs)
    τS = estimate_delay(s, opt.dmethod, τs)
    f = opt.method
    statT = f(t, τT, dims)
    statS = f(s, τS, dims)
    dT = dims[argmin(statT)]
    dS = dims[argmin(statS)]

    return EmbeddingTE(; dT = dT, τT = -τT, dS = dS, τS = -τS, ηTf = opt.η, dTf = 1)
end

function optimize_marginals(opt::OptimiseTraditional, s, t, c)
    τs = getlags(opt, t)
    dims = 1:opt.maxdim
    τT = estimate_delay(t, opt.dmethod, τs)
    τS = estimate_delay(s, opt.dmethod, τs)
    τC = estimate_delay(c, opt.dmethod, τs)

    f = opt.method
    statC = f(c, τC, dims)
    statT = f(t, τT, dims)
    statS = f(s, τS, dims)

    dC = dims[argmin(statC)]
    dT = dims[argmin(statT)]
    dS = dims[argmin(statS)]

    return EmbeddingTE(; dT = dT, τT = -τT, dS = dS, τS = -τS, dC = dC, τC = -τC,
        ηTf = opt.η, dTf = 1) # always predict a one-dimensional target vector.


end

EmbeddingTE(opt::OptimiseTraditional, x...) = optimize_marginals(opt, x...)
