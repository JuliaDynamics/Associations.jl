export RenyiCrossEntropy

"""
    RenyiCrossEntropy <: CrossEntropyEstimator

A nearest neighbor-based estimator of Rényi differential cross entropy that
uses a variant of LoftsGaarden & Quesenberry (1965)'s kNN-based density estimates
to compute probabilities. It doesn't do any bias correction.

!!! note "A novel estimator?"
    We've not been able to locate this estimator in the literature, so it might be new.
    Have you found it in a paper? Please let us know.
"""
Base.@kwdef struct RenyiCrossEntropy <: CrossEntropyEstimator
    k::Int = 1
    w::Int = 0
end

function entropy_cross(e::Renyi, est::RenyiCrossEntropy,
        x::AbstractDataset{D}, y::AbstractDataset{D}) where D
    e.q != 1 || error("`entropy_cross` not defined for `RenyiCrossEntropy` estimator for Renyi with q = $(e.q)")
    n, m = length(x), length(y)

    (; k, w) = est
    n, m = length(x), length(y)
    treex = KDTree(x, Euclidean())
    treey = KDTree(y, Euclidean())
    # Radii of the balls
    νs = last.(bulksearch(treey, x, NeighborNumber(k), Theiler(w))[2])
    bv = Entropies.ball_volume(D) # unit ball volume
    a = log.((bv .* νs .^ D .* m) ./ exp(digamma(k)))
    return a
    @show a .^ (e.q - 1)
    #c = 1/n * sum( a .^ (e.q - 1) )
    #c  = digamma(k) +
    #    log(bv) +
    #    log(m) + (D / n) * sum(log.(νs))
    #c = sum(log(k / (m*bv*ν^D)) for ν in νs)

    #c = sum((k / (m*bv*ν^D))^(e.q-1) for (ρ, ν) in zip(ρs, νs))
    #c /= n
    #return (1 / (1 - e.q) * c) / log(e.base, ℯ)
end

export RenyiCrossEntropyLord
Base.@kwdef struct RenyiCrossEntropyLord <: CrossEntropyEstimator
    k::Int = 1
    w::Int = 0
end

function entropy_cross(e::Renyi, est::RenyiCrossEntropyLord, x::AbstractDataset{D}, y::AbstractDataset{D}) where {D}
    # TODO: only for Shannon()
    (; k, w) = est
    N = length(x)
    treex = KDTree(x, Euclidean())
    treey = KDTree(y, Euclidean())
    idxsx, dsx = bulksearch(treex, y, NeighborNumber(k), Theiler(w))
    idxsy, dsy = bulksearch(treey, x, NeighborNumber(k), Theiler(w))

    # Decrease allocations and speed up computations by pre-allocating.
    # We're only dealing with matrices which in either axis has a maximum dimension of `k`,
    # so this is still far within the realm of where StaticArrays shines.
    # -------------------------------------------------------------------------------------
    # Contains neighborhood-centroid-centered vectors, where
    # `Cx[1]` := the centered query point
    # `Cx[1 + j]` := the centered `j`-th neighbor of the query point.
    Cx = MVector{k + 1, MVector{D}}(@SVector zeros(D) for i = 1:k+1)
    Cy = MVector{k + 1, MVector{D}}(@SVector zeros(D) for i = 1:k+1)

    # Centered neighbors need to be ordered row-wise in a matrix. We re-fill this matrix
    # for every query point `xᵢ`
    Ax = @MMatrix zeros(k+1, D)
    Ay = @MMatrix zeros(k+1, D)

    rsx = zeros(D)
    rsy = zeros(D)

    c = 0.0

    for (i, (xᵢ, yᵢ)) in enumerate(zip(x, y))
        neighborsᵢx = @views x[idxsx[i]]
        neighborsᵢy = @views y[idxsy[i]]

        # Center neighborhood around mean of the neighborhood.
        cx = centroid(xᵢ, neighborsᵢx, Cx)
        cy = centroid(yᵢ, neighborsᵢy, Cy)
        center_neighborhood!(cx, Cx, xᵢ, neighborsᵢx) # put centered vectors in `Cx`
        center_neighborhood!(cy, Cy, yᵢ, neighborsᵢy) # put centered vectors in `Cy`
        fill_A!(Ax, Cx, D)
        fill_A!(Ay, Cy, D)

        # SVD. The columns of Vt are the semi-axes of the ellipsoid, while Σx and Σx give
        # the magnitudes of the axes in each marginal.
        Ux, Σx, Vtx = svd(Ax)
        Uy, Σy, Vty = svd(Ay)

        σ₁x = Σx[1]
        σ₁y = Σy[1]

        ϵᵢx = last(dsx[i])
        ϵᵢy = last(dsy[i])

        # Scale semi-axis lengths to k-th neighbor distance
        rsx .= ϵᵢx .* (Σx ./ σ₁x)
        rsy .= ϵᵢy .* (Σy ./ σ₁y)

        # Create matrix ellipse representation, centered at origin.
        Λx = hyperellipsoid_matrix(Vtx, rsx)
        Λy = hyperellipsoid_matrix(Vty, rsy)

        # Take the point `xᵢ` as the origin for the neighborhood, so we can check directly
        # whether points are inside the ellipse. This happens for a point `x`
        # whenever `xᵀΛx <= 1`.
        nns_centeredx = (pt - xᵢ for pt in neighborsᵢx)
        nns_centeredy = (pt - yᵢ for pt in neighborsᵢy)

        kᵢx = count([transpose(px) * Λx * px <= 1.0 for px in nns_centeredx])
        kᵢy = count([transpose(py) * Λy * py <= 1.0 for py in nns_centeredy])

        # In the paper the point itself is always counted inside the ellipsoid,
        # so that there is always one point present. Here we instead set the local density
        # to zero (i.e. skip this iteration step) if that is the case.
        if kᵢx > 0 && kᵢy > 0
            volᵢx = gamma(1 + D/2) / (N * π^(D/2)) * prod(Σx ./ σ₁x) * ϵᵢx^D
            volᵢy = gamma(1 + D/2) / ((N - 1) * π^(D/2)) * prod(Σy ./ σ₁y) * ϵᵢy^D
            f̂xᵢ = (kᵢx / N) / volᵢx
            ĝyᵢ = (kᵢy / N) / volᵢy
            c += log(f̂xᵢ * ĝyᵢ^(e.q - 1))
        end
    end
    c /= N
    c *= -1 / (1 - e.q)
    return c / log(e.base, ℯ)
end

export RenyiDifferentialCrossEntropy
struct RenyiDifferentialCrossEntropy <: CrossEntropyDefinition end



# Analytical expressions for Thierrin et al. (2022)'s definition, which they call
# (discrete) *Rényi cross entropy*. Note: this definition is different from what
# they call the *natural Rényi cross entropy, for which analytical expressions appear
# in a separate file.
# --------------------------------------------------------------------
function entropy_cross(e::Renyi, definition::RenyiDifferentialCrossEntropy,
        x::Exponential, y::Exponential)
    @assert e.q != 1
    q, b = e.q, e.base
    λ₁, λ₂  = x.θ, y.θ
    λₕ = λ₁ + (e.q - 1)*λ₂
    c = 1 / (1 - e.q) * log(λ₁ / λₕ) - log( λ₂)
    return c / log(b, ℯ)
end

function entropy_cross(e::Renyi, definition::RenyiDifferentialCrossEntropy,
        x::Normal, y::Normal)
    q, b = e.q, e.base
    @assert e.q != 1
    σ₁, σ₂ = x.σ, y.σ
    μ₁, μ₂ = x.μ, y.μ
    σ²ₕ = σ₂^2 + (q - 1)*σ₁^2
    @assert σ²ₕ > 0

    c = 0.5*(log(2π * σ₂^2) + 1/(1 - q)*log(σ₂^2 / σ²ₕ) + ((μ₁ - μ₂)^2)/σ²ₕ)
    return c / log(b, ℯ)
end

function entropy_cross(e::Renyi, definition::RenyiDifferentialCrossEntropy,
    x::MvNormal, y::MvNormal)
    q, b = e.q, e.base
    @assert e.q != 1
    Σ₁, Σ₂ = x.Σ, y.Σ
    μ₁, μ₂ = x.μ, y.μ
    Nx, Ny = length(μ₁), length(μ₂)
    @assert Nx == Ny; N = Nx

    μ₁ᵀ, μ₂ᵀ = transpose(μ₁), transpose(μ₂)
    Σ₁⁻¹, Σ₂⁻² = inv(Σ₁), inv(Σ₂)
    A = Σ₁⁻¹ + (q - 1) * Σ₂⁻²
    A⁻¹ = inv(A)
    d = μ₁ᵀ * Σ₁⁻¹ * μ₁ +
        (q - 1) * μ₂ᵀ * Σ₂⁻² * μ₂ -
        (μ₁ᵀ * Σ₁⁻¹ + (q - 1) *  μ₂ᵀ * Σ₂⁻²) * A⁻¹ * (Σ₁⁻¹*μ₁ + (q - 1)*Σ₂⁻²*μ₂)
    c = 1 / (2 - 2*q) * (-log(det(A)*det(Σ₁)) + (1 - q)*log(2π)^N * det(Σ₂) - d)
    return c / log(b, ℯ)
end
