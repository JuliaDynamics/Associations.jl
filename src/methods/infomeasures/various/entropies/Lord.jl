using StaticArrays: @MVector, @MMatrix, SVector, MVector
using Neighborhood: NeighborNumber, KDTree, Theiler, Euclidean
using Neighborhood: bulksearch
using LinearAlgebra: svd

export Lord

"""
    Lord <: DifferentialEntropyEstimator
    Lord(k = 1, w = 0)

`Lord` estimates the [`Shannon`](@ref) differential entropy using a nearest neighbor
approach with a local nonuniformity correction (LNC).

## Description

Assume we have samples ``\\bar{X} = \\{\\bf{x}_1, \\bf{x}_2, \\ldots, \\bf{x}_N \\}`` from a
continuous random variable ``X \\in \\mathbb{R}^d`` with support ``\\mathcal{X}`` and
density function ``f : \\mathbb{R}^d \\to \\mathbb{R}``. `Lord` estimates the
[Shannon](@ref) differential entropy

```math
H(X) = \\int_{\\mathcal{X}} f(x) \\log f(x) dx = \\mathbb{E}[-\\log(f(X))],
```

by using the resubstitution formula

```math
\\hat{\\bar{X}, k} = -\\mathbb{E}[\\log(f(X))]
\\approx \\sum_{i = 1}^N \\log(\\hat{f}(\\bf{x}_i)),
```

where ``\\hat{f}(\\bf{x}_i)`` is an estimate of the density at ``\\bf{x}_i`` constructed
in a manner such that ``\\hat{f}(\\bf{x}_i) \\propto \\dfrac{k(x_i) / N}{V_i}``,
where ``k(x_i)`` is the number of points in the neighborhood of ``\\bf{x}_i``, and ``V_i``
is the volume of that neighborhood.

While most nearest-neighbor based differential entropy estimators uses regular volume
elements (e.g. hypercubes, hyperrectangles, hyperspheres) to approximate local densities,
the `Lord` estimator uses hyperellopsoid volume elements to estimate the local density.
These hyperellipsoids are, for each query point `xᵢ`, estimated using singular value
decomposition (SVD) on the `k`-th nearest neighbors of `xᵢ`. Thus, the hyperellipsoids
stretch/compress in response to the local geometry around each sample point. This
makes `Lord` a well-suited entropy estimator for a wide range of systems.

[^Lord2015]:
    Lord, W. M., Sun, J., & Bollt, E. M. (2018). Geometric k-nearest neighbor estimation
    of entropy and mutual information. Chaos: An Interdisciplinary Journal of Nonlinear
    Science, 28(3), 033114.
"""
Base.@kwdef struct Lord <: DifferentialEntropyEstimator
    k::Int = 10
    w::Int = 0
end

import ComplexityMeasures: entropy
function entropy(e::Shannon, est::Lord, x::AbstractDataset{D}) where {D}
    (; k, w) = est
    N = length(x)
    tree = KDTree(x, Euclidean())
    knn_idxs, ds = bulksearch(tree, x, NeighborNumber(k), Theiler(w))

    # Decrease allocations and speed up computations by pre-allocating.
    # We're only dealing with matrices which in either axis has a maximum dimension of `k`,
    # so this is still far within the realm of where StaticArrays shines.
    # -------------------------------------------------------------------------------------
    rs = @MVector zeros(D) # Scaled ellipsoid axis lengths
    Λ = @MMatrix zeros(D, D) # Hyperellipsoid matrix

    # C contains neighborhood-centroid-centered vectors, where
    # `C[1]` := the centered query point
    # `C[1 + j]` := the centered `j`-th neighbor of the query point.
    C = MVector{k + 1, MVector{D}}(@SVector zeros(D) for i = 1:k+1)
    # Centered neighbors need to be ordered row-wise in a matrix. We re-fill this matrix
    # for every query point `xᵢ`
    A = @MMatrix zeros(k+1, D)

    # Precompute some factors
    γ = gamma(1 + D/2)
    f = N * π^(D/2)

    h = 0.0
    for (i, xᵢ) in enumerate(x)
        neighborsᵢ = @views x[knn_idxs[i]]
        # Center neighborhood points around mean of the neighborhood.
        c = centroid(xᵢ, neighborsᵢ, C)
        center_neighborhood!(C, c, xᵢ, neighborsᵢ) # put centered vectors in `C`
        fill_A!(A, C, D)

        # SVD. The columns of Vt are the semi-axes of the ellipsoid, while Σ gives the
        # magnitudes of the axes.
        U, Σ, Vt = svd(A)
        σ₁ = Σ[1]
        ϵᵢ = last(ds[i])

        # Scale semi-axis lengths to k-th neighbor distance
        rs .= ϵᵢ .* (Σ ./ σ₁)
        # Create matrix ellipse representation, centered at origin (fill Λ)
        hyperellipsoid_matrix!(Λ, Vt, rs)
        # Take the point `xᵢ` as the origin for the neighborhood, so we can check directly
        # whether points are inside the ellipse. This happens for a point `x`
        # whenever `xᵀΛx <= 1`.
        nns_centered = (pt - xᵢ for pt in neighborsᵢ)
        kᵢ = count([transpose(p) * Λ * p <= 1.0 for p in nns_centered])
        # In the paper the point itself is always counted inside the ellipsoid,
        # so that there is always one point present. Here we instead set the local density
        # to zero (by just skipping the computation) if that is the case.
        if kᵢ > 0
            h += log(kᵢ * γ / (f * ϵᵢ^D * prod(Σ ./ σ₁)) )
        end
    end
    h = - h / N

    return h / log(e.base, ℯ)
end
entropy(est::Lord, args...) = entropy(Shannon(), est, args...)

function fill_A!(A, C, D)
    for (j, m) in enumerate(C)
        A[j, :] .= SVector{D}(m)
    end
end

"""
    center_neighborhood!(c, C, xᵢ, neighbors)

Center the point `xᵢ`, as well as each of its neighboring points `nⱼ ∈ neighbors`,
to the (precomputed) centroid `c` of the points `{xᵢ, n₁, n₂, …, nₖ}`, and store the
centered vectors in the pre-allocated vector `C`.
"""
function center_neighborhood!(C::MVector{K, V}, c, xᵢ,
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

function hyperellipsoid_matrix!(Λ, directions, extents)
    Λ .= 0.0
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
