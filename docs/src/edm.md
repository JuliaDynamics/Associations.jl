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

```@example simplex_projection
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

### Example: determining optimal embedding dimension

```@docs
delay_simplex
```

The simplex projection method can also be used to determine the optimal embedding dimension for a time series.
Given an embedding lag `τ`, we can embed a time series `x` for a range of embedding dimensions `d ∈ 2:dmax` and
compute the average prediction power over multiple `ks` using the simplex projection method.

Here, we compute the average prediction skills from `k=1` up to `k=15` time steps into the future, for 
embedding dimensions `d = 2:10`.

```@example
using CausalityTools, DynamicalSystems, Plots; gr()

sys = CausalityTools.ExampleSystems.lorenz_lorenz_bidir()
T, Δt = 150, 0.05
lorenz = trajectory(sys, T, Δt = Δt, Ttr = 100)[:, 1:3]
x1, x2, x3 = columns(lorenz)

# Determine the optimal embedding delay
τ = estimate_delay(x1, "ac_zero")

# Compute average prediction skill for 
ds, ks = 2:10, 1:15
ρs = delay_simplex(x1, τ, ds = ds, ks = ks)

plot(xlabel = "Embedding dimension", ylabel = "ρ̄(observed, predicted")
plot!(ds, ρs, label = "", c = :black, marker = :star)
savefig("simplex_embedding.svg") # hide
```

![](simplex_embedding.svg)

Based on the predictability criterion, the optimal embedding dimension, for this particular realization
of the first variable of the Lorenz system, seems to be 2.

## S-map

```@docs
smap
```

The s-map, or sequential locally weighted global map, was introduced in Sugihara (1994)[^Sugihara1994]. The s-map approximates the dynamics of a system as a locally weighted global map, with a tuning parameter ``\theta`` that controls the degree of nonlinearity in the model. For ``\theta = 0``, the model is the maximum likelihood global linear solution (of eq. 2 in Sugihara, 1994), and for increasing ``\theta > 0``, the model becomes increasingly nonlinear and localized (Sugihara, 1996)[^Sugihara1996].

When such a model has been constructed, it be used as prediction tool for out-of-sample points. It extends the prediction power of the simplex projection method, which only approximates the dynamics with a piecewise linear map. Let's elaborate with an example.

### Example: prediction power for the Lorenz system

In our implementation of `smap`, the input is a multivariate dataset - which can be a `Dataset` of either the raw variables of a multivariate dynamical system, or a `Dataset` containing an embedding of a single time series. 
Here, we'll show an example of the former.

Let's generate an example orbit from a bidirectionally coupled set of Lorenz systems. We'll use the built-in 
`CausalityTools.ExampleSystems.lorenz_lorenz_bidir` system, and select the first three variables for analysis.

```@example smap_lorenz
using CausalityTools, Plots, DynamicalSystems, Statistics; gr()

sys = CausalityTools.ExampleSystems.lorenz_lorenz_bidir()
T, Δt = 150, 0.05
lorenz = trajectory(sys, T, Δt = Δt, Ttr = 100)[:, 1:3]
p_orbit = plot(xlabel = "x", ylabel = "y", zlabel = "z")
plot!(columns(lorenz)..., marker = :circle, label = "", ms = 2, lα = 0.5)
savefig("lorenz_attractor.svg") # hide
```

![](lorenz_attractor.svg)

Now, we compute the `k`-step forward predictions for `k` ranging from `1` to `15`. The tuning parameter `θ` 
varies from `0.0` (linear model) to `2.0` (strongly nonlinear model). Our goal is to see which model 
yields the best predictions across multiple `k`.

We'll use the first 500 points of the orbit to train the model. Then, using that model, we try to
predict the next 500 points (which are not part of the training set). 
Finally, we compute the correlation between the predicted values and the observed values, which measures
the model prediction skill. This procedure is repeated for each combination of `k` and `θ`.

```@example smap_lorenz
ks, θs = 1:15, 0.0:0.5:2.0
n_trainees, n_predictees = 500, 500;

# Compute correlations between predicted values `preds` and actual values `truths` 
# for all parameter combinations
cors = zeros(length(ks), length(θs))
for (i, k) in enumerate(ks)
    for (j, θ) in enumerate(θs)
        pred = n_trainees+1:n_trainees+n_predictees
        train = 1:n_trainees

        preds, truths = smap(lorenz, θ = θ, k = k, trainees = train, predictees = pred)
        cors[i, j] = cor(preds, truths)
    end
end

p_θ_k_sensitivity = plot(xlabel = "Prediction time (k)", ylabel = "cor(observed, predicted)", 
    legend = :bottomleft, ylims = (-1.1, 1.1))
hline!([0], ls = :dash, c = :grey, label = "")
markers = [:star :circle :square :star5 :hexagon :circle :star]
cols = [:black, :red, :blue, :green, :purple, :grey, :black]
labels = ["θ = $θ" for θ in θs]
for i = 1:length(θs)
    plot!(ks, cors[:, i], marker = markers[i], c = cols[i], ms = 5, label = labels[i])
end

# Let's also plot the subset of points we're actually using.
p_orbit = plot(columns(lorenz[1:n_trainees+n_predictees])..., marker = :circle, label = "", ms = 2, lα = 0.5)

plot(p_orbit, p_θ_k_sensitivity, layout = grid(1, 2), size = (700, 400))
savefig("smap_sensitivity_lorenz.svg") # hide
```

![](smap_sensitivity_lorenz.svg)

The nonlinear models (colored lines and symbols) far outperform the linear model (black line + stars).

Because the predictions for our system improves with increasingly nonlinear models, it indicates that our 
system has some inherent nonlinearity. This is, of course, correct, since our Lorenz system is chaotic.

A formal way to test the presence of nonlinearity is, for example, to define the null hypothesis 
"H0: predictions do not improve when using an equivalent nonlinear versus a linear model" (equivalent in the sense  that the only parameter is `θ`) or, equivalently, "`H0: ρ_linear = ρ_nonlinear`". If predictions do in fact improve, we
instead accept the alternative hypothesis that prediction *do* improve when using nonlinear models versus using linear models. This can be formally tested using a z-statistic [^Sugihara1994].

[^Sugihara1994]: Sugihara, G. (1994). Nonlinear forecasting for the classification of natural time series. Philosophical Transactions of the Royal Society of London. Series A: Physical and Engineering Sciences, 348(1688), 477-495.
[^Sugihara1996]: Sugihara, George, et al. "Nonlinear control of heart rate variability in human infants." Proceedings of the National Academy of Sciences 93.6 (1996): 2608-2613.