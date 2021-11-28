# Empirical dynamical modelling

## Simplex projection

```@docs
simplex_predictions
```

### Example: reproducing Sugihara & May (1990)

The simplex projection method was introduced in Sugihara & May (1990). Let's try to reproduce their figure 1.

We start by defining the tent map and generating a time series for `μ = 1.98`.

```@example simplex_projection
using CausalityTools, DynamicalSystems, Plots, Statistics, Distributions; gr()

function eom_tentmap(dx, x, p, n)
    x = x[1]
    μ = p[1]
    dx[1] = x < 0.5 ? μ*x : μ*(1 - x)

    return
end

function tentmap(u₀ = rand(); μ = 1.98)
    DiscreteDynamicalSystem(eom_tentmap, [u₀], [μ])
end 

npts = 2000
sys = tentmap(μ = 1.98)
ts = diff(trajectory(sys, npts , Ttr = 1000)[:, 1])

p_ts = plot(xlabel = "Time step", ylabel = "Value", ylims = (-1.1, 1.1))
plot!(ts[1:100], label = "")
savefig("simplex_ts.svg") # hide
```

![](simplex_ts.svg)

Next, let's compute the predicted and observed values for `k = 1` and `k = 7` using embedding dimension 3 and 
embedding lag 1. We'll use the first 1000 points as our training set, and try to predict the next 500 points. 

```@example simplex_projection
d, τ = 3, 1
training = 1:1000
prediction = 1001:1500
x_1, x̃_1 = simplex_predictions(ts, 1, d = d, τ = τ, training = training, prediction = prediction)
x_7, x̃_7 = simplex_predictions(ts, 7, d = d, τ = τ, training = training, prediction = prediction)

p_obs_vs_pred = plot(xlabel = "Observed values", ylabel = "Predicted values")
scatter!(x_1, x̃_1, label = "k = 1", shape = :circle)
scatter!(x_7, x̃_7, label = "k = 7", shape = :star5)
savefig("simplex_correspondence.svg") # hide
```

![](simplex_correspondence.svg)

There is high correlation between observed and predicted values when predicting only one time step (`k = 1`)
into the future. As `k` increases, the performance drops off. Let's investigate this systematically.

```@example simplex_projection
kmax = 20
cors = zeros(kmax)
for k = 1:kmax
    X, X̃ = simplex_predictions(ts, k, d = d, τ = τ, training = training, prediction = prediction)
    cors[k] = cor(X, X̃)
end

plot(legend = :bottomleft, ylabel = "Correlation coefficient (ρ)", 
    xlabel = "Prediction time (k)",
    ylims = (-1.1, 1.1))
hline!([0], ls = :dash, label = "", color = :grey)
scatter!(1:kmax, cors, label = "")
savefig("simplex_correlation.svg") # hide
```

![](simplex_correlation.svg)

The correlation between observed and predicted values is near perfect until `k = 3`, and then rapidly 
drops off as `k` increases. At `k = 8`, there is virtually no correlation between observed and predicted values.
This means that, for this particular system, for this particular choice of embedding and choice of training/prediction sets, the predictability of the system is limited to about 4 or 5 time steps into the future (if you want good predictions). 

The main point of Sugihara & May's paper was that this drop-off of prediction accuracy with `k` is characteristic of chaotic systems, and can be used to distinguish chaos from regular behaviour in time series.

Let's demonstrate this by also investigating how the correlation between observed and predicted values behaves as a function of `k` for a regular, non-chaotic time series. We'll use a sine wave with additive noise.

```@example
𝒩 = Uniform(-0.5, 0.5)
xs = 0.0:1.0:2000.0
r = sin.(0.5 .* xs) .+ rand(𝒩, length(xs))
plot(r[1:200])

cors_sine = zeros(kmax)
for k = 1:kmax
    X, X̃ = simplex_predictions(r, k, d = d, τ = τ, training = training, prediction = prediction)
    cors_sine[k] = cor(X, X̃)
end

plot(legend = :bottomleft, ylabel = "Correlation coefficient (ρ)", 
    xlabel = "Prediction time (k)",
    ylims = (-1.1, 1.1))
hline!([0], ls = :dash, label = "", color = :grey)
scatter!(1:kmax, cors, label = "tent map", marker = :star)
scatter!(1:kmax, cors_sine, label = "sine")
savefig("simplex_correlation_sine_tent.svg") # hide
```

![](simplex_correlation_sine_tent.svg)

In contrast to the tent map, for which prediction accuracy drops off and stabilizes around zero for increasing `k`, the prediction accuracy is rather insensitive to the choice of `k` for the noisy sine time series. 