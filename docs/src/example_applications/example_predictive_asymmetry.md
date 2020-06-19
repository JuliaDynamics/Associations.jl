# Predictive asymmetry

Here, we demonstrate how the predictive asymmetry method on the tent map example from Schreiber (2000)'s 
transfer entropy paper.

## [`VisitationFrequency`](@ref) estimator

```@example
using CausalityTools, DynamicalSystems, Random, Plots, Statistics, LaTeXStrings; pyplot()
Random.seed!(1234);

function ulam(dx, x, p, t)
    f(x) = 2 - x^2
    ε = p[1]
    dx[1] = f(ε*x[length(dx)] + (1-ε)*x[1])
    for i in 2:length(dx)
        dx[i] = f(ε*x[i-1] + (1-ε)*x[i])
    end
end

# Initialize system and trigger compilation by generating a trajectory.
n = 500 # time series length
D = 30 # number of variables
ds = DiscreteDynamicalSystem(ulam, rand(D) .- 0.5, [0.04])
tr = trajectory(ds, n; Ttr = 1000);

# Use VisitationFrequency estimator
method = VisitationFrequency(binning = ExtendedPalusLimit())

# Define prediction lags (we're guessing that a maximum predction lag of 8 is good;
# thus ignoring any information transfer beyond this time horizon)
ηs = 1:8

εs = 0.0:0.05:1.0
pa_x1x2 = zeros(D-1, length(εs)); 
pa_x2x1 = zeros(D-1, length(εs))

for (k, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, n; Ttr = 1000)

    for i in 1:D-1
        X1 = tr[:, i]; X2 = tr[:, i+1]
        @assert !any(isnan, X1)
        @assert !any(isnan, X2)

        # Pick 𝒜(η = max(ηs)) to plot
        #@show k, i, X1[1:5], X2[1:5]
        pa_x1x2[i, k] = predictive_asymmetry(X1, X2, ηs, method)[end]
        pa_x2x1[i, k] = predictive_asymmetry(X2, X1, ηs, method)[end]
    end
end

# Plot the predictive asymmetries for each of the D-1 pairs of variables (x_i, x_{i+1})
plot(xlabel = L"\epsilon", ylabel = L"\mathcal{A}", legend = :right,
    fg_legend = :transparent, bg_legend = :transparent, size = (382, 300),
    guidefont = font(11), tickfont = font(11), legendfont = font(11))

plot!(εs, pa_x1x2[1, :], label = L"x_i \to x_{i+1}", c = :black, lα = 0.5, lw = 0.5)
plot!(εs, pa_x2x1[1, :], label = L"x_{i+1} \to x_i", c = :red, lα = 0.5, lw = 0.5)

for i = 2:D-1
    plot!(εs, pa_x1x2[i, :], label = "", c = :black, lα = 0.5, lw = 0.5)
    plot!(εs, pa_x2x1[i, :], label = "", c = :red, lα = 0.5, lw = 0.5)
end
hline!([1], ls = :dash, label = "", c = :grey)
```


## [`TransferOperatorGrid`](@ref) estimator

```@example
using CausalityTools, DynamicalSystems, Random, Plots, Statistics, LaTeXStrings; pyplot()
Random.seed!(1234);

function ulam(dx, x, p, t)
    f(x) = 2 - x^2
    ε = p[1]
    dx[1] = f(ε*x[length(dx)] + (1-ε)*x[1])
    for i in 2:length(dx)
        dx[i] = f(ε*x[i-1] + (1-ε)*x[i])
    end
end

# Initialize system and trigger compilation by generating a trajectory.
n = 500 # time series length
D = 30 # number of variables
ds = DiscreteDynamicalSystem(ulam, rand(D) .- 0.5, [0.04])
tr = trajectory(ds, n; Ttr = 1000);

# Use TranferOperatorGrid estimator
method = TransferOperatorGrid(binning = ExtendedPalusLimit())

# Define prediction lags (we're guessing that a maximum predction lag of 8 is good;
# thus ignoring any information transfer beyond this time horizon)
ηs = 1:5

εs = 0.0:0.05:1.0
pa_x1x2 = zeros(D-1, length(εs)); 
pa_x2x1 = zeros(D-1, length(εs))

for (k, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, n; Ttr = 1000)

    for i in 1:D-1
        X1 = tr[:, i]; X2 = tr[:, i+1]
        @assert !any(isnan, X1)
        @assert !any(isnan, X2)

        # Pick 𝒜(η = max(ηs)) to plot
        #@show k, i, X1[1:5], X2[1:5]
        pa_x1x2[i, k] = predictive_asymmetry(X1, X2, ηs, method)[end]
        pa_x2x1[i, k] = predictive_asymmetry(X2, X1, ηs, method)[end]
    end
end

# Plot the predictive asymmetries for each of the D-1 pairs of variables (x_i, x_{i+1})
plot(xlabel = L"\epsilon", ylabel = L"\mathcal{A}", legend = :right,
    fg_legend = :transparent, bg_legend = :transparent, size = (382, 300),
    guidefont = font(11), tickfont = font(11), legendfont = font(11))

plot!(εs, pa_x1x2[1, :], label = L"x_i \to x_{i+1}", c = :black, lα = 0.5, lw = 0.5)
plot!(εs, pa_x2x1[1, :], label = L"x_{i+1} \to x_i", c = :red, lα = 0.5, lw = 0.5)

for i = 2:D-1
    plot!(εs, pa_x1x2[i, :], label = "", c = :black, lα = 0.5, lw = 0.5)
    plot!(εs, pa_x2x1[i, :], label = "", c = :red, lα = 0.5, lw = 0.5)
end
hline!([1], ls = :dash, label = "", c = :grey)
```

## [`SymbolicPerm`](@ref) estimator

```@example
using CausalityTools, DynamicalSystems, Random, Plots, Statistics, LaTeXStrings; pyplot()
Random.seed!(1234);

function ulam(dx, x, p, t)
    f(x) = 2 - x^2
    ε = p[1]
    dx[1] = f(ε*x[length(dx)] + (1-ε)*x[1])
    for i in 2:length(dx)
        dx[i] = f(ε*x[i-1] + (1-ε)*x[i])
    end
end

# Initialize system and trigger compilation by generating a trajectory.
n = 500 # time series length
D = 30 # number of variables
ds = DiscreteDynamicalSystem(ulam, rand(D) .- 0.5, [0.04])
tr = trajectory(ds, n; Ttr = 1000);

# Use motifs of length 2.
method = SymbolicPerm(b = 2, m = 2)

# Use default embedding settings (embedding delays -m+1, -m+2, ..., 0 , and prediction lags 1, ..., m
# for symbolic estimators)
E = EmbeddingTE()

# Define prediction lags (we're guessing that a maximum predction lag of 3*m=6 is good;
# thus ignoring any information transfer beyond this time horizon)
ηs = 1:3

εs = 0.0:0.05:1.0
pa_x1x2 = zeros(D-1, length(εs)); 
pa_x2x1 = zeros(D-1, length(εs))

for (k, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, n; Ttr = 1000)

    for i in 1:D-1
        X1 = tr[:, i]; X2 = tr[:, i+1]
        @assert !any(isnan, X1)
        @assert !any(isnan, X2)

        # Pick 𝒜(η = max(ηs)) to plot
        #@show k, i, X1[1:5], X2[1:5]
        pa_x1x2[i, k] = predictive_asymmetry(X1, X2, ηs, method)[end]
        pa_x2x1[i, k] = predictive_asymmetry(X2, X1, ηs, method)[end]
    end
end

# Plot the predictive asymmetries for each of the D-1 pairs of variables (x_i, x_{i+1})
plot(xlabel = L"\epsilon", ylabel = L"\mathcal{A}", legend = :right,
    fg_legend = :transparent, bg_legend = :transparent, size = (382, 300),
    guidefont = font(11), tickfont = font(11), legendfont = font(11))

plot!(εs, pa_x1x2[1, :], label = L"x_i \to x_{i+1}", c = :black, lα = 0.5, lw = 0.5)
plot!(εs, pa_x2x1[1, :], label = L"x_{i+1} \to x_i", c = :red, lα = 0.5, lw = 0.5)

for i = 2:D-1
    plot!(εs, pa_x1x2[i, :], label = "", c = :black, lα = 0.5, lw = 0.5)
    plot!(εs, pa_x2x1[i, :], label = "", c = :red, lα = 0.5, lw = 0.5)
end
hline!([1], ls = :dash, label = "", c = :grey)
```
