using Entropies: KozachenkoLeonenko, ball_volume
using Neighborhood
using SpecialFunctions

export crossentropy

# Find the distances to the k-th nearest neighbors.
function maximum_neighbor_distances_cross(x, y; w::Int = 0, k::Int = 2)
    theiler = Theiler(w)
    w >= 0 || error("w, the number of neighbors to exclude, must be >= 0")
    tree_y = KDTree(y, Euclidean())
    idxs, dists = bulksearch(tree_y, x, NeighborNumber(k), theiler)
    # Distances to k-th nearest neighbor of each of xᵢ ∈ x
    return [d[end] for d in dists]
end

"""
    crossentropy(x, y, est::KozachenkoLeonenko)

Compute the cross-entropy

```math
H_y(x) = \\sum{x_i} p(x) \\log{\\dfrac{1}{q(x)}},
```

of the `D`-dimensional [`Dataset`](@ref)s `x` and `y`.
"""
function crossentropy(x::AbstractDataset{D, T}, y::AbstractDataset{D, T},
        est::KozachenkoLeonenko) where {D, T}
    (; k, w, base) = eest
    Nx = length(x)
    Ny = length(y)
    xds = maximum_neighbor_distances_cross(x, y, k = k, w = est.w)
    bv = ball_volume(D) # for *Euclidean* metric.
    sum_term = sum(log.(est.k ./ (Ny * bv * (xds .^ D)) ) )
    return log(Nx - 1) + log(bv) - (digamma(est.k) + (D / Nx) * sum_term)
end

function bentropy(x::AbstractDataset{D, T}, est::KozachenkoLeonenko) where {D, T}
    (; k, w, base) = est
    Nx = length(x)
    xds = Entropies.maximum_neighbor_distances(x, w, k)
    bv = ball_volume(D) # for *Euclidean* metric.
    return - digamma(k) + log(Nx - 1) + log(base, bv) +
        ((D / Nx) * sum(log.(base, xds .^ D)))
end
