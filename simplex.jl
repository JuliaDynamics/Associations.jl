using CausalityTools, DynamicalSystems, Plots, StatsBase, Statistics, Distributions, Neighborhood, LinearAlgebra; gr()

function eom_tentmap(dx, x, p, n)
    x = x[1]
    μ = p[1]
    dx[1] = x < 0.5 ? μ*x : μ*(1 - x)

    return
end

function tentmap(u₀ = rand(); μ = 1.9)
    DiscreteDynamicalSystem(eom_tentmap, [u₀], [μ])
end 

npts = 2000
sys = tentmap(μ = 1.9)
ts = trajectory(sys, npts , Ttr = 1000)[:, 1]

plot(legend = :topleft)
k = 3
kmax = 20
cors = zeros(kmax)
τ = 1
d = 3
@show ts
training = 1:200
prediction = 201:500
for k = 1:kmax
    X, X̃ = simplex_predictions(ts, k, d = d, τ = τ, training = training, prediction = prediction)
    cors[k] = cor(X, X̃)
end

plot(legend = :bottomleft, ylabel = "Correlation coefficient (ρ)", 
    xlabel = "Prediction time (k)",
    ylims = (-0.05, 1.1))
scatter!(1:kmax, cors, label = "")



# r = 1:35
# p1 = plot(prediction[r], ts[prediction][r], c = :black, lw = 2, alpha = 0.5)
# #plot!(prediction_set, x, c = :green, lw = 2)
# plot!(prediction[r], xpred[r], ls = :dot, c = :blue)
# @show xpred

# r = 1:35
# p1 = plot(prediction, ts[prediction], c = :black, lw = 2, alpha = 0.5)
# #plot!(prediction_set, x, c = :green, lw = 2)
# plot!(prediction, xpred, ls = :dot, c = :blue)
