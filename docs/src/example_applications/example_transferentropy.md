# Transfer entropy

Here, we reproduce the tent map example from Schreiber (2000)[^Schreiber2000] using different 
transfer entropy estimators.

[^Schreiber2000]: Schreiber, Thomas. "Measuring information transfer." Physical review letters 85.2 (2000): 461.

## `VisitationFrequency` estimator 

```@example 
using TransferEntropy, DynamicalSystems, Random, Plots, LaTeXStrings; pyplot()
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
n = 300 # time series length
D = 100
ds = DiscreteDynamicalSystem(ulam, rand(D) .- 0.5, [0.04])
tr = trajectory(ds, n; Ttr = 1000);

# Compute transfer entropy between variables x1 and x2 for different values
# of the parameter ε
estimator = VisitationFrequency(b = 2)
E = EmbeddingTE()

εs = 0.0:0.01:1.0
te_x1x2 = zero(εs); te_x2x1 = zero(εs)
for (i, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, n; Ttr = 1000)
    X1 = tr[:, 1]; X2 = tr[:, 2]
    @assert !any(isnan, X1)
    @assert !any(isnan, X2)
    te_x1x2[i] = transferentropy(X1, X2, E, estimator)
    te_x2x1[i] = transferentropy(X2, X1, E, estimator)
end

plot(xlabel = L"\epsilon", ylabel = "TE (bits)", legend = :bottomright)
plot!(εs, te_x1x2, label = L"x_1 \to x_2")
plot!(εs, te_x2x1, ls = :dash, label = L"x_2 \to x_1")
```

## `TransferOperatorGrid` estimator 

```@example 
using TransferEntropy, DynamicalSystems, Random, Plots, LaTeXStrings; pyplot()
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
n = 300 # time series length
D = 100
ds = DiscreteDynamicalSystem(ulam, rand(D) .- 0.5, [0.04])
tr = trajectory(ds, n; Ttr = 1000);

# Compute transfer entropy between variables x1 and x2 for different values
# of the parameter ε
estimator = TransferOperatorGrid(b = 2)
E = EmbeddingTE()

εs = 0.0:0.01:1.0
te_x1x2 = zero(εs); te_x2x1 = zero(εs)
for (i, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, n; Ttr = 1000)
    X1 = tr[:, 1]; X2 = tr[:, 2]
    @assert !any(isnan, X1)
    @assert !any(isnan, X2)
    te_x1x2[i] = transferentropy(X1, X2, E, estimator)
    te_x2x1[i] = transferentropy(X2, X1, E, estimator)
end

plot(xlabel = L"\epsilon", ylabel = "TE (bits)", legend = :bottomright)
plot!(εs, te_x1x2, label = L"x_1 \to x_2")
plot!(εs, te_x2x1, ls = :dash, label = L"x_2 \to x_1")1)")
```

## `NearestNeighborMI` estimator 


```@example 
using TransferEntropy, DynamicalSystems, Random, Plots, LaTeXStrings; pyplot()
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
n = 300 # time series length
D = 100
ds = DiscreteDynamicalSystem(ulam, rand(D) .- 0.5, [0.04])
tr = trajectory(ds, n; Ttr = 1000);

# Compute transfer entropy between variables x1 and x2 for different values
# of the parameter ε
estimator = NearestNeighborMI(b = 2)
E = EmbeddingTE()

εs = 0.0:0.01:1.0
te_x1x2 = zero(εs); te_x2x1 = zero(εs)
for (i, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, n; Ttr = 1000)
    X1 = tr[:, 1]; X2 = tr[:, 2]
    @assert !any(isnan, X1)
    @assert !any(isnan, X2)
    te_x1x2[i] = transferentropy(X1, X2, E, estimator)
    te_x2x1[i] = transferentropy(X2, X1, E, estimator)
end

plot(xlabel = L"\epsilon", ylabel = "TE (bits)", legend = :bottomright)
plot!(εs, te_x1x2, label = L"x_1 \to x_2")
plot!(εs, te_x2x1, ls = :dash, label = L"x_2 \to x_1")
```


## [`SymbolicPerm`](@ref) estimator

The [`SymbolicPerm`](@ref) estimator is highly sensitive to the motif length `m`.
For the Ulam map example, the time series mostly jump up/down with a some intermittent
sections where values decrease or increase. We will therefore use the shortest 
possible motif length `m=2`.

```@example 
using TransferEntropy, DynamicalSystems, Random, Plots, LaTeXStrings; pyplot()
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
n = 300 # time series length
D = 100
ds = DiscreteDynamicalSystem(ulam, rand(D) .- 0.5, [0.04])
tr = trajectory(ds, n; Ttr = 1000);

# Use motifs of length 2.
estimator = SymbolicPerm(b = 2, m = 2)

# Use default embedding settings (embedding delays -m+1, -m+2, ..., 0 , and prediction lags 1, ..., m)
E = EmbeddingTE()

εs = 0.0:0.01:1.0
te_x1x2 = zero(εs); te_x2x1 = zero(εs)
for (i, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, n; Ttr = 1000)
    X1 = tr[:, 1]; X2 = tr[:, 2]
    @assert !any(isnan, X1)
    @assert !any(isnan, X2)
    te_x1x2[i] = transferentropy(X1, X2, E, estimator)
    te_x2x1[i] = transferentropy(X2, X1, E, estimator)
end

plot(xlabel = L"\epsilon", ylabel = "TE (bits)", legend = :bottomright)
plot!(εs, te_x1x2, label = L"x_1 \to x_2")
plot!(εs, te_x2x1, ls = :dash, label = L"x_2 \to x_1")
```


## `SymbolicAmplitudeAware` estimator 

The [`SymbolicAmplitudeAware`](@ref) estimator is also highly sensitive to the motif length `m`.
For the Ulam map example, the time series mostly jump up/down with a some intermittent
sections where values decrease or increase. We will therefore use the shortest 
possible motif length `m=2`.

```@example 
using TransferEntropy, DynamicalSystems, Random, Plots, LaTeXStrings; pyplot()
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
n = 300 # time series length
D = 100
ds = DiscreteDynamicalSystem(ulam, rand(D) .- 0.5, [0.04])
tr = trajectory(ds, n; Ttr = 1000);

# Use motifs of length 2 and A = 0.5, which weights amplitudes and 
# differences in consecutive amplitudes within marginal delay embedding vectors 
# equally.
estimator = SymbolicAmplitudeAware(b = 2, m = 2, A = 0.5)

# Use default embedding settings (embedding delays -m+1, -m+2, ..., 0 , and prediction lags 1, ..., m)
E = EmbeddingTE()

εs = 0.0:0.01:1.0
te_x1x2 = zero(εs); te_x2x1 = zero(εs)
for (i, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, n; Ttr = 1000)
    X1 = tr[:, 1]; X2 = tr[:, 2]
    @assert !any(isnan, X1)
    @assert !any(isnan, X2)
    te_x1x2[i] = transferentropy(X1, X2, E, estimator)
    te_x2x1[i] = transferentropy(X2, X1, E, estimator)
end

plot(xlabel = L"\epsilon", ylabel = "TE (bits)", legend = :bottomright)
plot!(εs, te_x1x2, label = L"x_1 \to x_2")
plot!(εs, te_x2x1, ls = :dash, label = L"x_2 \to x_1")
```

