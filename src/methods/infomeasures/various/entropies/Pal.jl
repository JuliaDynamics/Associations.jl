export Pal

"""
    Pal <: <: DifferentialEntropyEstimator
    Pal(; k, w, p)

The Pal [`Renyi`](@ref) differential entropy estimator.

[^Pál2010]:
    Pál, D., Póczos, B., & Szepesvári, C. (2010). Estimation of Rényi entropy and mutual
    information based on generalized nearest-neighbor graphs. Advances in Neural
    Information Processing Systems, 23.
"""
Base.@kwdef struct Pal{P} <: DifferentialEntropyEstimator
    k::Int = 1
    w::Int = 0
    p::P = 1
end

# V := A set of finite points in R^d
# NNs(V) := generalized nearest neighbor graph on V
# Edge set of NNs(V) contains for each i ∈ S a vertex from x ∈ V to its i-th nearest neighbor
# S := Finite set of positive integers
# k := maximum element of S

function Lₚ(est::Pal, v, p, x::AbstractDataset{D}) where D
    (; k, w, p) = est
    N = length(x)
    tree = KDTree(x, Euclidean())
    idxs, ds = bulksearch(tree, x, NeighborNumber(k), Theiler(w))
    Lₚ = 0.0
    for i = 1:N
        Lₚ += sum(ds[i])
    end
    return Lₚ
end
