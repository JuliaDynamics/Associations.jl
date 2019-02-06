## Different ways of partitioning
The `transferentropy_transferoperator_grid` and `transferentropy_freq`
estimators both operate on partitions on the state space. Below, we demonstrate
the four different ways of discretizing the state space.

```@setup TE_partitioning_schemes
using CausalityTools
using Plots
using Statistics
using TimeseriesSurrogates
```

First, let's create some example time series, embed them and organize the
computation of marginal probabilities.

```@example TE_partitioning_schemes
x = cumsum(rand(300))
y = sin.(cumsum(rand(300)))*0.3 .+ x

τ = 1 # embedding lag
ν = 1 # forward prediction lag
E_xtoy = customembed([x, y], [2, 2, 2, 1], [ν, 0, -τ, 0])
E_ytox = customembed([y, x], [2, 2, 2, 1], [ν, 0, -τ, 0])

# Organize marginals
Tf = [1]     # target, future
Tpp = [2, 3] # target, present and past
Spp = [4]    # source, present (and past, if we wanted)
v = TEVars(Tf, Tpp, Spp)
```


### Hyper-rectangles by subdivision of axes (`ϵ::Int`)

First, we use an integer number of subdivisions along each axis of the delay
embedding when partitioning.

```@example TE_partitioning_schemes
ϵs = 1:2:50 # integer number of subdivisions along each axis of the embedding
te_estimates_xtoy = zeros(length(ϵs))
te_estimates_ytox = zeros(length(ϵs))
v = TEVars([1], [2, 3], [4])

for (i, ϵ) in enumerate(ϵs)
    te_estimates_xtoy[i] = tetogrid(E_xtoy, ϵ, v)
    te_estimates_ytox[i] = tetogrid(E_ytox, ϵ, v)
end

plot(ϵs, te_estimates_xtoy, label = "TE(x -> y)", lc = :black)
plot!(ϵs, te_estimates_ytox, label = "TE(y -> x)", lc = :red)
xlabel!("# subdivisions along each axis")
ylabel!("Transfer entropy (bits)")
```

### Hyper-cubes of fixed size (`ϵ::Float`)
We do precisely the same, but use fixed-width hyper-cube bins. The values of the
logistic map take values on `[0, 1]`, so using bins width edge lengths `0.1`
should give a covering corresponding to using `10` subdivisions along
each axis of the delay embedding. We let `ϵ` take values on `[0.05, 0.5]`.

```@example TE_partitioning_schemes
ϵs = 0.02:0.02:0.5
te_estimates_xtoy = zeros(length(ϵs))
te_estimates_ytox = zeros(length(ϵs))
v = TEVars([1], [2, 3], [4])

for (i, ϵ) in enumerate(ϵs)
    te_estimates_xtoy[i] = tetogrid(E_xtoy, ϵ, v)
    te_estimates_ytox[i] = tetogrid(E_ytox, ϵ, v)
end

plot(ϵs, te_estimates_xtoy, label = "TE(x -> y)", lc = :black)
plot!(ϵs, te_estimates_ytox, label = "TE(y -> x)", lc = :red)
xlabel!("Hypercube edge length")
ylabel!("Transfer entropy (bits)")
xflip!()
```

### Hyper-rectangles of fixed size (`ϵ::Vector{Float}`)

It is also possible to use hyper-rectangles, by specifying the edge lengths
along each coordinate axis of the delay embedding. In our case, we use a
four-dimensional, embedding, so we must provide a 4-element vector of edge
lengths

```@example TE_partitioning_schemes
# Define slightly different edge lengths along each axis
ϵs_x1 = LinRange(0.05, 0.5, 10)
ϵs_x2 = LinRange(0.02, 0.4, 10)
ϵs_x3 = LinRange(0.08, 0.6, 10)
ϵs_x4 = LinRange(0.10, 0.3, 10)

te_estimates_xtoy = zeros(length(ϵs_x1))
te_estimates_ytox = zeros(length(ϵs_x1))
v = TEVars([1], [2, 3], [4])
mean_ϵs = zeros(10)
for i ∈ 1:10
    ϵ = [ϵs_x1[i], ϵs_x2[i], ϵs_x3[i], ϵs_x4[i]]
    te_estimates_xtoy[i] = tetogrid(E_xtoy, ϵ, v)
    te_estimates_ytox[i] = tetogrid(E_ytox, ϵ, v)

    # Store average edge length (for plotting)
    mean_ϵs[i] = mean(ϵ)
end

plot(mean_ϵs, te_estimates_xtoy, label = "TE(x -> y)", lc = :black)
plot!(mean_ϵs, te_estimates_ytox, label = "TE(y -> x)", lc = :red)
xlabel!("Average hypercube edge length")
ylabel!("Transfer entropy (bits)")
xflip!()
```


### Hyper-rectangles by variable-width subdivision of axes (`ϵ::Vector{Int}`)

Another way to construct hyper-rectangles is to subdivide each
coordinate axis into segments of equal length. In our case, we use a
four-dimensional, embedding, so we must provide a 4-element vector providing
the number of subdivisions we want along each axis.

```@example TE_partitioning_schemes
# Define different number of subdivisions along each axis.
ϵs = 3:50
mean_ϵs = zeros(length(ϵs))

te_estimates_xtoy = zeros(length(ϵs))
te_estimates_ytox = zeros(length(ϵs))
v = TEVars([1], [2, 3], [4])
for (i, ϵᵢ) ∈ enumerate(ϵs)
    ϵ = [ϵᵢ - 1, ϵᵢ, ϵᵢ, ϵᵢ + 1]
    te_estimates_xtoy[i] = tetogrid(E_xtoy, ϵ, v)
    te_estimates_ytox[i] = tetogrid(E_ytox, ϵ, v)

    # Store average number of subdivisions for plotting
    mean_ϵs[i] = mean(ϵ)
end

plot(mean_ϵs, te_estimates_xtoy, label = "TE(x -> y)", lc = :black)
plot!(mean_ϵs, te_estimates_ytox, label = "TE(y -> x)", lc = :red)
xlabel!("Average number of subdivisions along the embedding axes")
ylabel!("Transfer entropy (bits)")
```
