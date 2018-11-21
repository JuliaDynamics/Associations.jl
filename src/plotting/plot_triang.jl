using Simplices: even_sampling_rules
using StaticArrays
using RecipesBase

function interiorpoints(s::SArray{D, T}, n::Int) where {D, T}
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
getsimplex(r::Embeddings.AbstractEmbedding, t::DelaunayTriangulation, i::Int) =
    r.points[:, t.indices[:, i]]

getsimplex(d::Dataset, t::DelaunayTriangulation, i::Int) =
    r.points[:, t.indices[:, i]]
"""
    forwardmap(r::Embeddings.AbstractEmbedding,
                t::DelaunayTriangulation,
                i::Int)

Get the vertices of the i-th simplex of the triangulation projected
one step forward in time.
"""
forwardmap(r::Embeddings.AbstractEmbedding, t::DelaunayTriangulation, i::Int) =
    r.points[:, t.indices[:, i] .+ 1]

forwardmap(r::Dataset, t::DelaunayTriangulation, i::Int) =
    r.points[:, t.indices[:, i] .+ 1]

"""
    plot_triang(r::AbstractReconstruction, t::DelaunayTriangulation)

Plot the triangulation of a 3D state space reconstruction.
"""
function plot_triang(r::Embeddings.AbstractEmbedding, t::DelaunayTriangulation;
        plot_states = true,
        plot_simplices = true,
        evenly_subsample = true,
        k = 3, # The splitting factor, k^dim new points are introduced per simplex
        n_interior_pts = 100,
        insert_pts = false,
        ms = 0.5, mc = :red,
        edges = true,
        vertices = true,
        highlight = 0,
        highlight_ms = 2,
        extrapts = []) # highlight a given simplex, 0 = no simplices are highlighted

    n_simplices = size(t.indices, 2)
    n_points = size(r.points, 2)
    p = plot()

    # Embedding center
    #ce = sum(r.points, 2)/size(r.points, 2)
    #lp = r.points[:, end]

    #scatter3d!(p, splitaxes(ce), ms = 3)
    #scatter3d!(p, splitaxes(lp), ms = 3, mc = :black)

    if length(extrapts) > 0
        scatter3d!(p, splitaxes(extrapts), ms = 1.4, mc = :black, mα = 0.5)
    end

    for i = 1:n_points
        pᵢx = r.points[1, i]
        pᵢy = r.points[2, i]
        pᵢz = r.points[3, i]
        #scatter3d!(p, pᵢx, pᵢy, pᵢz,
        #    mc = :black, mα = 0.5, ms = 1)
    end


    for i = 1:n_simplices
        # Get the ith simplex and connect its vertices
        sᵢ  = getsimplex(r, t, i)
        cᵢ = connectvertices(sᵢ)

        if edges
            plot3d!(p, preparesimplex(cᵢ), lc = :black, lw = 0.5, lα = 0.5, legend = false)
        end

        if vertices
            scatter3d!(p, preparesimplex(sᵢ), mc = mc, ms = ms, mα = 0.5, legend = false)
        end

        if insert_pts
            if !evenly_subsample
                random_interior_pts = interiorpoints(sᵢ, n_interior_pts)
                scatter3d!(p, splitaxes(random_interior_pts),
                    mc = :red, ms = ms, lα = 0.5)
            elseif evenly_subsample
                scatter3d!(p, splitaxes(sᵢ * even_sampling_rules(3, k)),
                        mc = :red, ms = ms, lα = 0.5)
            end
        end
    end
    # Starting simplex
    if length(r) - 1 > highlight> 0
        S_orig = getsimplex(r, t, highlight)
        S_image = forwardmap(r, t, highlight)

        if edges
            plot3d!(preparesimplex(S_orig), label = "simplex $highlight", lw = 2, lc = :green)
            plot3d!(preparesimplex(S_image), label = "forward projection", lw = 2, lc = :red)
        end

        if vertices
            scatter3d!(preparesimplex(S_orig), label = "simplex $highlight", ms = highlight_ms, mc = :green)
            scatter3d!(preparesimplex(S_image), label = "forward projection", ms = highlight_ms, mc = :red)
        end
    end
    p
end

export plot_triang, interiorpoints, splitaxes, preparesimplex, forwardmap
