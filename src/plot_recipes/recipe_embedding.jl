@recipe function f(E::StateSpaceReconstruction.AbstractEmbedding;
        vars = [1, 2, 3],
        mc = :black, ms = 1, lc = :black, lw = 1)
    if size(E.points, 1) < 3
        error("Cannot plot reconstruction with < 3 variables.")
    end
    legend --> false
    xlabel --> "x1"
    ylabel --> "x2"
    zlabel --> "x3"
    markercolor --> mc
    markersize --> ms
    ([E.points[i, :] for i in vars]...,)
end
