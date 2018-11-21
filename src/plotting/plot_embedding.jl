####################################
# Plotting
####################################
@recipe function f(r::AbstractEmbedding; dims = (1, 2, 3))
    pts = r.points

    if dimension(r) > 3
        @warn "Embedding dim > 3, plotting three first axes"
        pts = r.points[1:3, :]
    end

    if length(dims) > 3
        error("dim = $dim. Can be at most 3.")
    end

    if length(dims) == 3
        X = pts[dims[1], :]
        Y = pts[dims[2], :]
        Z = pts[dims[3], :]
        return X, Y, Z
    elseif length(dims) == 2
        X = pts[dims[1], :]
        Y = pts[dims[2], :]
        return X, Y
    end
end
