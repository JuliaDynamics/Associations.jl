#export PCPA

"""
    PCPA <: GraphAlgorithm
    PCPA(; est::E = FPVP(), τS::TS = Pecuzal(), τC::TC = Pecuzal(), ηT::ET = 1:5,
        α::A = 0.05, maxdepth::M = 1)

An algorithm for inferring causal graps from time series, based on the [`PA`](@ref) measure.

Used with [`infer_graph`](@ref) to infer a directed graph between the input variables.
This graph does not encode information about causal lags; it only indicates for which
variables there is evidence of dynamical coupling, given the time-asymmetric
predictability of the time series.
"""
Base.@kwdef struct PCPA{E, TS, TC, ET, A, M} <: GraphAlgorithm
    est::E = FPVP()
    τS::TS = Pecuzal()
    τC::TC = Pecuzal()
    ηT::ET = 1:5
    α::A = 0.05
    maxdepth::M = 1
end

function infer_graph(alg::PCPA, x)
    (; est, τS, τC, ηT, α, maxdepth) = alg
    m = PA(; τS, τC, ηT)
    test = PATest(m, est)
    alg = PCRobust(test, test; α, maxdepth)
    return infer_graph(alg, x)
end
