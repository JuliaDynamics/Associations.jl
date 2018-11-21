
#########################################################
#########################################################
# Visualizing different rectangular binnings
#########################################################
#########################################################

#########################################################
# Helper functions to plot 3D rectangles.
#########################################################
"""
    rectangle3dpts(x, y, z, ϵx, ϵy, ϵz)

Compute the coordinates of the vertices of a 3D rectangle, given the
coordinates of the origin (`x`, `y`, and `z`) and the corresponding edge lengths
(`ϵx`, `ϵy`, `ϵz`).
"""
function rectangle3dpts(x, y, z, ϵx, ϵy, ϵz)
    v1b = [x,    y,    z]
    v2b = [x+ϵx, y,    z]
    v3b = [x+ϵx, y+ϵy, z]
    v4b = [x,    y+ϵy, z]
    v1t = [x,    y,    z+ϵz]
    v2t = [x+ϵx, y,    z+ϵz]
    v3t = [x+ϵx, y+ϵy, z+ϵz]
    v4t = [x,    y+ϵy, z+ϵz]
    (v1b, v2b, v3b, v4b, v1t, v2t, v3t, v4t)
end

"""
    rectangle3dpts(o, ϵ)

Given the origin `o` (3-element vector) and edge lengths `ϵ` (also a
3-element vector) of a 3D rectangle, return the coordinates of its
vertices.
"""
function rectangle3dpts(o, ϵ)
    x, y, z = (o...,)
    ϵx, ϵy, ϵz = (ϵ...,)
    rectangle3dpts(x, y, z, ϵx, ϵy, ϵz)
end

"""
    connectvertices(o, ϵ)

Given the origin `o` (3-element vector) and edge lengths `ϵ` (also a
3-element vector) of a 3D rectangle, return a vector of line segments
(each a dim-by-n_vertices array) that when plotted together yields
the entire rectangle.
"""
function connectvertices(o, ϵ)
    # Get the vertices
    v1b, v2b, v3b, v4b, v1t, v2t, v3t, v4t = rectangle3dpts(o, ϵ)

    # Connect vertices in top and bottom planes
    connect_top = zeros(3, 5)
    connect_bottom = zeros(3, 5)

    # Connect corners of the top and bottom planes
    corner1 = zeros(3, 2)
    corner2 = zeros(3, 2)
    corner3 = zeros(3, 2)
    corner4 = zeros(3, 2)
    connect_bottom[:, 1:5] = [v1b v2b v3b v4b v1b]
    connect_top[:, 1:5]    = [v1t v2t v3t v4t v1t]
    corner1[:, 1:2] = [v1b v1t]
    corner2[:, 1:2] = [v2b v2t]
    corner3[:, 1:2] = [v3b v3t]
    corner4[:, 1:2] = [v4b v4t]

    return [connect_top, connect_bottom,
            corner1, corner2, corner3, corner4]
end

"""
    splitaxes(x)

Return a vector of the individual components of an array of points
provided as an array where each point is a column.
"""
splitaxes(x) = ([x[k, :] for k = 1:size(x, 1)]...,)

"""
    plot_3D_rect!(p, origin, edgelengths;
             lc = :black, lw = 1.0, ls = :solid)

Append a a 3D rectangle to an existing plot `p`, given the
origin `o` (3-element vector) and edge lengths
`ϵ` (3-element vector) of the rectangle.
"""
function plot_3D_rect!(p, o, ϵ;
             lc = :black, lw = 1.0, ls = :solid, lα = 0.7)

    linesegments = connectvertices(o, ϵ)
    for segment in linesegments
        plot!(p, splitaxes(segment),
            lc = lc, lw = lw, ls = ls, lα = lα)
    end
end


#########################################################
# Given a partition and a binning scheme given by ϵ,
# plot the grid superimposed on the points of the orbit.
#########################################################
"""
    plot_partition(pts::AbstractArray{T, 2}, ϵ;
                    mc = :blue, ms = 1.5, mα = 0.8,
                    lc = :black, lw = 2, ls = :dash, lα = 0.6) where T

Partition the space defined by `pts` into rectangular boxes
with a binning scheme controlled by `ϵ`.

The following `ϵ` will work:

* `ϵ::Int` divide each axis into `ϵ` intervals of the same size.
* `ϵ::Float` divide each axis into intervals of size `ϵ`.
* `ϵ::Vector{Int}` divide the i-th axis into `ϵᵢ` intervals of the same size.
* `ϵ::Vector{Float64}` divide the i-th axis into intervals of size `ϵᵢ`.

The points are assumed to be provided as an array where each point is
a column.

`mc`, `ms` and `mα` control the marker color, marker size and marker opacity,
respectively. `lc`, `lw`, `lα` and `ls` control the line color, line width,
line opacity and line style, respectively.
"""
function plot_partition(pts::AbstractArray{T, 2}, ϵ;
                mc = :black, ms = 2, mα = 0.8,
                lc = :blue, lw = 1.5, ls = :dash, lα = 0.6) where T

    if size(pts, 1) > size(pts, 2)
        #info("Treating each row as a point")
        pts = transpose(pts)
    end

    axisminima, stepsizes = minima_and_stepsizes(pts, ϵ)

    # Assign bins to the points.
    v = assign_coordinate_labels(pts, ϵ)

    # Bins may be visited by several times, so we'll get rid
    # of the repetitions.
    V = unique(v, dims = 2)
    n_visited_boxes = size(V, 2)

    # Plot the partition grid over the points of the reconstructed
    # orbit.
    p = plot(legend = false)
    for i = 1:n_visited_boxes
        origin = V[:, i]
        plot_3D_rect!(p, origin, stepsizes, lc = lc, ls = ls, lα = lα)
    end
    scatter3d!(p, splitaxes(pts), ms = ms, mc = mc, mα = mα)
    xlabel!("x")
    ylabel!("y")
    p
end

"""
    plot_partition(E::Embeddings.AbstractEmbedding, ϵ; vars = [1, 2, 3],
                    mc = :blue, ms = 2, mα = 0.8,
                    lc = :black, lw = 2, ls = :dash, lα = 0.6)

Partition the embedding into rectangular boxes with a binning scheme controlled
by `ϵ`. If there are more than three variables in the embedding, you can set
which one to use with the `vars` argument (by default, `vars = [1, 2, 3]`).

The following `ϵ` will work:

* `ϵ::Int` divide each axis into `ϵ` intervals of the same size.
* `ϵ::Float` divide each axis into intervals of size `ϵ`.
* `ϵ::Vector{Int}` divide the i-th axis into `ϵᵢ` intervals of the same size.
* `ϵ::Vector{Float64}` divide the i-th axis into intervals of size `ϵᵢ`.


`mc`, `ms` and `mα` control the marker color, marker size and marker opacity,
respectively. `lc`, `lw`, `lα` and `ls` control the line color, line width,
line opacity and line style, respectively.
"""
function plot_partition(E::Embeddings.AbstractEmbedding, ϵ;
                vars = [1, 2, 3],
                mc = :blue, ms = 2, mα = 0.8,
                lc = :black, lw = 1.5, ls = :dash, lα = 0.6)
    plot_partition(E.points[vars, :], ϵ, mc = mc, ms = ms,
                    lc = lc, lw = lw, ls = ls)
end

export plot_partition, plot_3D_rect!, splitaxes, connectvertices
