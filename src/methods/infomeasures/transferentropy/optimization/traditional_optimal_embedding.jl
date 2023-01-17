using DelayEmbeddings

export OptimiseTraditional
export optimize_marginals_te

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
    optimize_marginals_te([scheme = OptimiseTraditional()], s, t, [c])

Optimize marginal embeddings for transfer entropy computation from source time series `s`
to target time series `t`, conditioned on `c` if `c` is given, using the provided
optimization `scheme`.
"""
function optimize_marginals_te end
getlags(opt::OptimiseTraditional{<:Float64}, x) = 1:floor(Int, length(x)*opt.maxlag)
getlags(opt::OptimiseTraditional{<:Int}, args...) = 1:opt.maxlag

"""
    optimize_marginals_te(opt::OptimiseTraditional, s, t, [c]; exclude_source = false) → EmbeddingTE

Optimise the marginals for a transfer entropy analysis from source time series `s` to
target time series `t`, potentially given a conditional time series `c`.

If `exclude_source == true`, then no optimisation is done for the source. This is
useful for [`SurrogateTest`](@ref), because most surrogate methods accept
univariate time series, and if we embed the source and it becomes multidimensional,
then we can't create surrogates. A future optimization is to do column-wise surrogate
generation.
"""
function optimize_marginals_te(opt::OptimiseTraditional, s, t; exclude_source = false)
    τs = getlags(opt, t)
    dims = 1:opt.maxdim
    f = opt.method

    if exclude_source
        τT = estimate_delay(t, opt.dmethod, τs)
        statT = f(t, τT, dims)
        dT = dims[argmin(statT)]
        return EmbeddingTE(; dT = dT, τT = -τT, ηTf = opt.η, dTf = 1)
    else
        τT = estimate_delay(t, opt.dmethod, τs)
        τS = estimate_delay(s, opt.dmethod, τs)
        statT = f(t, τT, dims)
        statS = f(s, τS, dims)
        dT = dims[argmin(statT)]
        dS = dims[argmin(statS)]
        return EmbeddingTE(; dT = dT, τT = -τT, dS = dS, τS = -τS, ηTf = opt.η, dTf = 1)
    end
end

function optimize_marginals_te(opt::OptimiseTraditional, s, t, c; exclude_source = false)
    τs = getlags(opt, t)
    dims = 1:opt.maxdim
    if exclude_source
        τT = estimate_delay(t, opt.dmethod, τs)
        τC = estimate_delay(c, opt.dmethod, τs)
        f = opt.method
        statC = f(c, τC, dims)
        statT = f(t, τT, dims)
        dC = dims[argmin(statC)]
        dT = dims[argmin(statT)]
        return EmbeddingTE(; dT = dT, τT = -τT, dC = dC, τC = -τC, ηTf = opt.η, dTf = 1) # always predict a one-dimensional target vector.
    else
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
end

EmbeddingTE(opt::OptimiseTraditional, x...) = optimize_marginals_te(opt, x...)
