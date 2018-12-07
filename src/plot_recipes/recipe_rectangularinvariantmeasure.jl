@recipe function f(riv::PerronFrobenius.RectangularInvariantMeasure;
                    boxfillfactor::Int = 3, linesegments = true,
                    lw = 0.8, lc = :black, lα = 0.5, ls = :do,
                    ms = 0.2, mc = :black, mα = 0.3)

    pts = riv.pts
    ϵ = riv.ϵ
    axisminima, stepsizes = minima_and_stepsizes(pts, ϵ)

    # All bins are assigned a bin origin.
    v = unique(riv.visited_bins_coordinates, dims = 2)

    # Bins may be visited by several times, so we'll get rid
    # of the repetitions.
    n_visited_boxes = size(v, 2)
    D = size(v, 1)

    # Plot the partition grid over the points of the reconstructed
    # orbit.
    legend --> false

    for i = 1:n_visited_boxes
        origin = riv.visited_bins_coordinates[:, i]
        # Get the measure of the box, and fill the box with a
        # number of points scaled to the measure.
        μ = riv.measure.dist[i]
        fillpts = fill_hyperrectangle_3D(origin,
                    stepsizes,
                    ceil(Int, μ*100)*boxfillfactor)
        @series begin
            seriestype := :scatter3d
            markercolor --> mc
            markeralpha --> μ*5
            markersize --> μ*10
            linecolor --> mc
            fillcolor --> mc
            fillalpha --> μ*2
            splitaxes(fillpts)
        end

        if linesegments
            linesegments = connectvertices(origin, stepsizes)

            for segment in linesegments
                @series begin
                    seriestype := :path
                    linealpha --> lα
                    linewidth --> lw
                    linecolor --> lc
                    splitaxes(segment)
                end
            end
        end
    end
end
