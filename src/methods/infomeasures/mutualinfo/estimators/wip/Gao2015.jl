
"""
    Gao2015 <: MutualInformationEstimator
    Gao2015(k = 1, w = 0, base = 2)

`Gao2015` estimates the [`Shannon`](@ref) mutual information using a nearest neighbor
approach with a local nonuniformity correction (LNC).

[Gao2015]
    Gao, S., Ver Steeg, G. &amp; Galstyan, A.. (2015). Efficient Estimation of Mutual
    Information for Strongly Dependent Variables. *Proceedings of the Eighteenth
        International Conference on Artificial Intelligence and Statistics*, in
        *Proceedings of Machine Learning Research* 38:277-286.
        Available from https://proceedings.mlr.press/v38/gao15.html.
"""
Base.@kwdef struct Gao2015{M} <: MutualInformationEstimator
    k::Int = 1
    w::Int = 0
    metric::M = Euclidean()
end

function mutualinfo(e::Renyi, est::Gao2015{B, Chebyshev}, x) where {B}
    (; k, w) = est

    joint = Dataset(x...)
    N = length(joint)
    M = length(x)

    # For each point `xᵢ ∈ joint`, find the index of its k-th nearest neighbor.
    tree_joint = KDTree(joint, metric_joint)
    idxs_joint = bulkisearch(tree_joint, joint, NeighborNumber(k + 1), Theiler(w))


end


function lnc_correction(est::Gao2015, x::AbstractDataset{D}, idxs_neighbors)
    (; k, w, base) = est

end

using Statistics
using LinearAlgebra

function pca(xᵢ, neighbors::SubDataset{D}) where D
    μ = xᵢ # manually set mean to be xᵢ, so that it is at the center of rotated rectangle.
    M = Matrix(x)
    C = @SMatrix cov(x)
    E = eigen(C)
    E.vectors[:, sortperm(E.values)]

end
