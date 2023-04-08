export PCMCIPlus

Base.@kwdef struct PCMCIPlus{P, C}
    τmax::Int = 1
    α::Real = 0.05
    pairwise_test::P = CorrTest()
    conditional_test::C = CorrTest()
end

Base.@kwdef mutable struct LaggedParents{PJ, PT}
    i::Int
    js::PJ = Vector{Int}(undef, 0)
    τs::PT = Vector{Int}(undef, 0)
end

function Base.show(io::IO, x::LaggedParents)
    s = ["x$(x.parents_js[i])($(x.parents_τs[i]))" for i in eachindex(x.parents)]
    all = "x$(x.i)(0) ← $(join(s, ", "))"
    show(io, all)
end

function infer_graph(alg::PCMCIPlus, X::AbstractStateSpaceSet)
    return infer_graph(alg, columns(X))
end

function infer_graph(alg::PCMCIPlus, X)
    parents = skeleton_lagged(alg, X)
    sc = skeleton_contemporaneous(alg, X, parents)
    return sc
end

include("skeleton_lagged.jl")
include("skeleton_contemp.jl")
