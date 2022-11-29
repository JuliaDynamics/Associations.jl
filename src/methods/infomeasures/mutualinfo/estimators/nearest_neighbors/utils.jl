using Distances: evaluate


# In the Kraskov1 estimator, ϵs are the distances in the Z = (X, Y) joint space
# In the Kraskov2 estimator, ϵs are the distances in the X and Y marginal spaces
function count_within_radius!(p, x, metric, ϵs, N)
    @inbounds for i in 1:N
        ϵ = ϵs[i] / 2
        xᵢ = x[i]
        p[i] = count(evaluate(metric, xᵢ, x[j]) < ϵ for j in 1:N)
    end

    return p
end
