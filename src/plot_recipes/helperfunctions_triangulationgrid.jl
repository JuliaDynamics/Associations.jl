using Simplices: subsample_coeffs
using Simplices: simplexintersection
using Simplices: intersectingvertices

function interiorpoints(s::SArray{D, T}, n::Int) where {D, T}
    R = rand(size(s, 2), n)
    convex_coeffs = ((1 ./ sum(R, dims = 1)) .* R)
    s * convex_coeffs
end

function interiorpoints(s, n::Int)
    R = rand(size(s, 2), n)
    convex_coeffs = ((1 ./ sum(R, dims = 1)) .* R)
    s * convex_coeffs
end

connectvertices(s) = hcat(s, hcat([s[:, i]
                    for i in [1, 3, 1, 4, 2]]...,))
splitaxes(x) = ([x[k, :] for k = 1:size(x, 1)]...,)
preparesimplex(s) = splitaxes(connectvertices(s))

"""
    getsimplex(r::Embeddings.AbstractEmbedding,
                t::DelaunayTriangulation,
                i::Int)

Get the vertices of the i-th simplex of the triangulation.
"""
getsimplex(E::Embeddings.AbstractEmbedding, DT::DelaunayTriangulation, i::Int) =
    E.points[:, DT.indices[:, i]]

getsimplex(d::Dataset, DT::DelaunayTriangulation, i::Int) =
    d.points[:, DT.indices[:, i]]

getsimplex(pts, DT::DelaunayTriangulation, i::Int) =
    pts[:, DT.indices[:, i]]

"""
    forwardmap(r::Embeddings.AbstractEmbedding,
                t::DelaunayTriangulation,
                i::Int)

Get the vertices of the i-th simplex of the triangulation projected
one step forward in time.
"""
forwardmap(r::Embeddings.AbstractEmbedding, t::DelaunayTriangulation, i::Int) =
    E.points[:, DT.indices[:, i] .+ 1]

forwardmap(r::Dataset, DT::DelaunayTriangulation, i::Int) =
    E.points[:, DT.indices[:, i] .+ 1]

forwardmap(pts, DT::DelaunayTriangulation, i::Int) =
    pts[:, DT.indices[:, i] .+ 1]



"""
Generate a total of `n_fillpoints` points lying inside the intersection
`sᵢ ∩ sⱼ` between simplices `sᵢ` and `sⱼ`.

This is done by first performing a Delaunay triangulation of the set of vertices
forming the intersecting volume. Then, points are inserted into the subsimplices
of the convex hull of the intersection.

If `sample_randomly = true`, then points are generated randomly inside each
subsimplex. If `sample_randomly = false` (default), then points are generated
as the centroids of the subsimplices resulting from a shape-preserving
refinement of the simplices in the convex hull.

## Arguments
- **``sᵢ``**: An array where each column is a simplex vertex.
- **``sⱼ``**: An array where each column is a simplex vertex.
- **``n_fillpoints``**: The total number of points to fill the convex
    hull of the intersection with. If the intersection is more complex than
    a single simplex, then the convex hull of the intersection is triangulated
    and `n_fillpoints` points are distributed evenly among the resulting
    subsimplices.
- **``sample_randomly``**: Should points be sampled randomly within the
    intersection? Default is `sample_randomly = false`, which inserts points
    inside the intersecting volume in a regular manner.
"""
function generate_fillpoints(sᵢ, sⱼ, n_fillpoints = 100, sample_randomly = true)

    intvol = simplexintersection(sᵢ, sⱼ)

    if intvol > 0 + 1e-9
        intersecting_vertices = intersectingvertices(sᵢ, sⱼ)
        if size(intersecting_vertices, 1) > 4
            triang = hcat(delaunay(intersecting_vertices)...,)
        else
            triang = zeros(Int, 4, 1)
            triang[:, 1] = 1:4
        end

        n_simplices = size(triang, 2)
        if n_simplices == 1
            npts_per_simplex = n_fillpoints
        elseif n_simplices > 1

            npts_per_simplex = ceil(Int, n_fillpoints/n_simplices)
        end
        C = transpose(subsample_coeffs(3, npts_per_simplex, sample_randomly))

        fillpts = Vector{Array{Float64, 2}}(undef, n_simplices)
        for k = 1:n_simplices
            fillpts[k] = transpose(C * intersecting_vertices[triang[:, k], :])
        end
        fillpts = hcat(fillpts...,)
        return intersecting_vertices, fillpts
    else
        return nothing
    end
end

export
interiorpoints,
connectvertices,
splitaxes,
preparesimplex,
getsimplex,
forwardmap,
generate_fillpoints
