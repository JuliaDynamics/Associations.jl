using StaticArrays: @MVector, @MMatrix, SVector, MVector
using Neighborhood: NeighborNumber, KDTree, Theiler, Euclidean
using Neighborhood: bulksearch

export Lord

"""
    Lord <: EntropyEstimator
    Lord(k = 1, w = 0, base = 2)

`Lord` estimates the [`Shannon`](@ref) mutual information using a nearest neighbor
approach with a local nonuniformity correction (LNC).

"""
Base.@kwdef struct Lord{M} <: EntropyEstimator
    k::Int = 10
    w::Int = 0
    metric::M = Euclidean()
end

import Entropies: entropy
function entropy(e::Renyi, est::Lord, x::AbstractDataset{D}) where {D}
    # TODO: only for Shannon()
    (; k, w, metric) = est
    N = length(x)
    tree = KDTree(x, metric)
    knn_idxs, ds = bulksearch(tree, x, NeighborNumber(k), Theiler(w))

    # Decrease allocations and speed up computations by pre-allocating.
    # We're only dealing with matrices which in either axis has a maximum dimension of `k`,
    # so this is still far within the realm of where StaticArrays shines.
    # -------------------------------------------------------------------------------------
    # Contains neighborhood-centroid-centered vectors, where
    # `C[1]` := the centered query point
    # `C[1 + j]` := the centered `j`-th neighbor of the query point.
    C = MVector{k + 1, MVector{D}}(@SVector zeros(D) for i = 1:k+1)

    # Centered neighbors need to be ordered row-wise in a matrix. We re-fill this matrix
    # for every query point `xᵢ`
    A = @MMatrix zeros(k+1, D)

    h = 0.0
    rs = zeros(D)
    g2 = gamma(D/2 + 1)
    πs = π^(D/2)
    ks = zeros(N)
    ϵs = zeros(N)
    ratios = zeros(N)
    for (i, xᵢ) in enumerate(x)
        neighborsᵢ = @views x[knn_idxs[i]]
        # Center neighborhood around mean of the neighborhood.
        c = centroid(xᵢ, neighborsᵢ, C)
        center_neighborhood!(c, C, xᵢ, neighborsᵢ) # put centered vectors in `M`
        fill_A!(A, C, D)

        # SVD. The columns of Vt are the semi-axes of the ellipsoid, while Σ gives the
        # magnitudes of the axes.
        U, Σ, Vt = svd(A)

        # How many of the `k` neighbors of `xᵢ` are within the ellipsoid
        # whose semi-axes have lengths rs[1], rs[2], ..., rs[D]?
        # ---------------------------------------------------------------------
        σ₁ = Σ[1]
        ϵᵢ = last(ds[i])
        # Scale semi-axis lengths to k-th neighbor distance, and create
        # matrix ellipse representation, relative to origin. Center points
        # around `xᵢ` too, so we can use the ellipse directly.
        rs .= ϵᵢ .* (Σ ./ σ₁)
        Λ = hyperellipsoid_matrix(Vt, rs)
        nns_centered = (pt - xᵢ for pt in neighborsᵢ)
        # points inside the ellipse obey xᵀΛx <= 1
        ks[i] = count([transpose(p) * Λ * p <= 1.0 for p in nns_centered]) + 1
        ratios[i] = sum(log.(Σ ./ σ₁))
        ϵs[i] = ϵᵢ
        h += log(ks[i] * gamma(1 + D/2) / (N * π^(D/2) * ϵᵢ^D * prod(Σ ./ σ₁)) )
    end
    h = - h / N
    # h = log(N) +
    #     log(π^(D / 2) / gamma(D/2 + 1)) -
    #     (1/N)*sum(log.(ks)) +
    #     (1/N)*sum(log.(ϵs) +
    #     (1/N)*sum(ratios) # already log-ed inside loop

    return h / log(e.base, ℯ)
end
entropy(est::Lord, args...; base = 2) = entropy(Shannon(; base), est, args...)

function fill_A!(A, C, D)
    for (j, m) in enumerate(C)
        A[j, :] .= SVector{D}(m)
    end
end

"""
    center_neighborhood!(c, C, xᵢ, neighbors)

Center the point `xᵢ`, as well as each of its neighboring points `nⱼ ∈ neighbors`,
to the (precomputed) centroid `c` of the points `{xᵢ, n₁, n₂, …, nₖ}`, and store the centered vectors
in the pre-allocated vector `C`.
"""
function center_neighborhood!(c, C::MVector{K, V}, xᵢ,
        neighbors::AbstractDataset{D}) where {V, K, D}
    rezero!(C)
    C[1] .= xᵢ .- c
    for (i, nᵢ) in enumerate(neighbors)
        C[1 + i] .= nᵢ .- c
    end
    return C
end

function rezero!(C)
    for mᵢ in C
        mᵢ .= 0.0
    end
end

function hyperellipsoid_matrix(directions, extents, D = length(extents))
    Λ = @MMatrix zeros(D, D)
    for i in axes(directions, 2)
        vᵢ = directions[:, i]
        Λ .+= (vᵢ * transpose(vᵢ)) ./ extents[i]^2
    end
    return Λ
end

function centroid(xᵢ, neighbors::AbstractDataset{D}, C) where D
    L = length(C) + 1
    centroid = MVector{D}(xᵢ)
    for nᵢ in neighbors
        centroid .+= nᵢ
    end
    centroid ./= L
    return SVector{D}(centroid)
end
