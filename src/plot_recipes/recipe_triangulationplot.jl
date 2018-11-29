include("helperfunctions_triangulationgrid.jl")

@recipe function plot_triang(original_pts, E::Embeddings.AbstractEmbedding, DT::DelaunayTriangulation;
            plot_states = true,
            plot_simplices = true,
            evenly_subsample = true,
            k = 3, # The splitting factor, k^dim new points are introduced per simplex
            n_interior_pts = 100,
            insert_pts = false,
            ms = 0.5,
            ms_originalpts = 1.5,
            mc = :black,
            mc_projection = :red,
            mc_originalpts = :black,
            mα = 0.5,
            edges = true,
            vertices = false,
            si = 0, # index of the i'th highlighted simplex
            sj = 0, # index of the j'th highlighted simplex
            ind_projection = 0, #show the simplex corresponding to the forward projection of the `ind_projection`-th simplex
            intersecting_vertices = nothing,
            plot_intersecting_pts = false,
            highlight_ms = 2,
            extrapts = [],
            lc = :black,
            plot_originalpts = true,
            plot_triang = true,
            plot_imagetriang = true,
            plot_fulltriang = false,
            lc_fulltriang = :black,
            lw_fulltriang = 0.5,
            lα_fulltriang = 0.5,
            lc_triang = :black,
            lc_imagetriang = :black,
            lw_imagetriang = 0.5,
            lw_triang = 0.5,
            ind_excluded_pt = 0,
            lw_highlight = 1.2,
            markershape_excluded_pt = :star5,
            mc_excluded_pt = :black,
            ms_excluded_pt = 3,
            lα = 0.2,
            lw = 0.5,
            lc_i = :green,
            lc_j = :red,
            lw_i = 0.5,
            lw_j = 0.5,
            lα_i = 1.0,
            lα_j = 1.0,
            mc_j = :red,
            mc_i = :green,
            mw_i = 0.5,
            mw_j = 0.5,
            mα_i = 1.0,
            mα_j = 1.0,
            lc_projection = :blue,
            lw_projection = 2,
            lα_projection = 1,
            idx_intersect_simplex = si, # the index of the simplex for which to compute the intersection between the simplex and its forward image
            mc_intersecting_pts = :black,
            markershape_intersecting_pts = :pentagon,
            ms_intersecting_pts = 1,
            fill_intersecting_volume = true,
            n_fillpoints = 500,
            ms_fillpoints = 1.0,
            mα_fillpoints = 1.0,
            mc_fillpoints = :black,
            sample_randomly = true
    )
    fulltriang = delaunay(StateSpaceReconstruction.embed(original_pts))
    n_simplices_full = size(fulltriang.indices, 2)
    n_simplices = size(DT.indices, 2)
    n_points = size(E.points, 2)

    # Embedding center
    #ce = sum(r.points, 2)/size(r.points, 2)
    #lp = r.points[:, end]

    #scatter3d!(p, splitaxes(ce), ms = 3)
    #scatter3d!(p, splitaxes(lp), ms = 3, mc = :black)

    legend --> false
    xaxis --> false
    if length(extrapts) > 0
        @series begin
            seriestype := :scatter
            legend --> false
            mc --> mc
            ms --> ms
            mα --> mα
            splitaxes(extrapts)
        end
    end
    if plot_originalpts
        @series begin
            seriestype := :scatter
            mc --> mc_originalpts
            ms --> ms_originalpts
            splitaxes(original_pts)
        end
        if ind_excluded_pt > 0
            @series begin
                seriestype := :scatter
                mc --> mc_excluded_pt
                markershape --> markershape_excluded_pt
                ms --> ms_excluded_pt
                splitaxes(original_pts[:, ind_excluded_pt])
            end
        end
    end

    if plot_intersecting_pts && (ind_projection > 0)
        sᵢ = getsimplex(original_pts, DT, si)
        sⱼ = forwardmap(original_pts, DT, ind_projection)
        intvol = simplexintersection(sᵢ, sⱼ)

        if intvol > 0 + 1e-8
            intersecting_vertices, fillpoints = generate_fillpoints(sᵢ, sⱼ,
                            n_fillpoints, sample_randomly)

            @series begin
                seriestype := :scatter
                ms --> ms_intersecting_pts
                markershape --> markershape_intersecting_pts
                mc --> mc_intersecting_pts
                splitaxes(transpose(intersecting_vertices))
            end
            @series begin
                seriestype := :scatter
                ms --> ms_fillpoints
                mα --> mα_fillpoints
                mc --> mc_fillpoints
                lc --> mc_fillpoints
                fc --> mc_fillpoints
                splitaxes(fillpoints)
            end
        end
    end

    if plot_fulltriang
        for i = 1:n_simplices
            sᵢ = getsimplex(original_pts, fulltriang, i)

            @series begin
                seriestype := :path
                lc --> lc_fulltriang
                lw --> lw_fulltriang
                lα --> lα_fulltriang
                legend --> false
                preparesimplex(sᵢ)
            end
        end
    end

    for i = 1:n_simplices
        sᵢ = getsimplex(original_pts, DT, i)

        if edges
            if plot_triang
                @series begin
                    seriestype := :path
                    lc --> lc_triang
                    lw --> lw_triang
                    lα --> lα
                    legend --> false
                    preparesimplex(sᵢ)
                end
            end

            if plot_imagetriang
                @series begin
                    seriestype := :path
                    lc --> lc_imagetriang
                    lw --> lw_imagetriang
                    lα --> lα
                    ls --> :dash
                    legend --> false
                    preparesimplex(forwardmap(original_pts, DT, i))
                end
            end
        end

        if vertices
            if plot_triang
                @series begin
                    seriestype := :scatter
                    mc --> mc
                    ms --> ms
                    mα --> mα
                    legend --> false
                end
            end
        end

        if insert_pts
            if !evenly_subsample
                @series begin
                    mc --> mc
                    ms --> ms
                    lα --> lα
                    legend --> false
                    random_interior_pts = interiorpoints(sᵢ, n_interior_pts)
                    splitaxes(random_interior_pts)
                end
            elseif evenly_subsample
                @series begin
                    mc --> mc
                    ms --> ms
                    lα --> lα
                    legend --> false
                    splitaxes(sᵢ * even_sampling_rules(3, k))
                end
            end
        end
    end

    if si > 0
        sᵢ = getsimplex(original_pts, DT, si)

        if edges
            @series begin
                seriestype := :path
                label --> "simplex #$sᵢ"
                lw --> lw_highlight
                lc --> lc_i
                lα --> lα_i
                preparesimplex(sᵢ)
            end
        end

        if vertices
            @series begin
                label --> "simplex #$sᵢ"
                ms --> ms_i
                mc --> mc_i
                preparesimplex(sᵢ)
            end
        end
    end

    if sj > 0
        sⱼ = getsimplex(original_pts, DT, sj)

        if edges
            @series begin
                seriestype := :path
                label --> "simplex #$sⱼ"
                lw --> lw_highlight
                lc --> lc_j
                lα --> lα_j
                ls --> :dot
                preparesimplex(sⱼ)
            end
        end

        if vertices
            @series begin
                label --> "simplex #$sⱼ"
                ms --> ms_j
                mc --> mc_j
                preparesimplex(sⱼ)
            end
        end
    end

    if ind_projection > 0
        @series begin
            seriestype := :path
            label --> "projection of simplex #$ind_projection"
            lw --> lw_highlight
            lc --> lc_projection
            lα --> lα_projection
            preparesimplex(forwardmap(original_pts, DT, ind_projection))
        end
    end
end


"""
Visualize a set of N points (`original_pts`), a delaunay triangulation
(`DT`) constructed from and embedding (`E`) of the first ``N-1`` points
of the points. By assuming a triangulation of the ``N-1`` first points,
the recipe allows for plotting both the original triangulation and the
image of that triangulation under the forward linear map of its vertices,
one step forward in time.

## Arguments
- **``original_pts``**: A set of ``N`` points where each column is a point.
- **``E``**: An embedding constructed from the first ``N-1`` points of
    `original_pts`.
- **``DT``**: A DelaunayTriangulation instance constructed from `E`.

## Keyword arguments
- **``plot_triang``**: Plot the original triangulation?
- **``plot_imagetriang``**: Plot the image of the triangulation under the forward linear
    map of its vertices?
- **``plot_fulltriang``**: Plot a triangulation of all points in `original_pts`?

"""
plot_triang

export plot_triang
