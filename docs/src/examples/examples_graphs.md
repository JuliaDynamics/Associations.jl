# Causal graphs

The causal graph [API](@ref causal_graphs) relies on the
[Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) package. Each
[`GraphAlgorithm`](@ref) returns some graph that is compatible with Graphs.jl.
To 

## [`PCRobust`](@ref)

The PC algorithm is perhaps the most famous algorithm for inferring causal graphs.
Here, we demonstrate the [`PCRobust`](@ref) variant on some random (uncoupled)
variables.

```@example causalgraph_corr
using CausalityTools
using Random
using Graphs
using CairoMakie, GraphMakie

# A function that plots an `n`-variable directed graph.
function plotgraph(g; labels)
    with_theme(theme_minimal(), resolution = (300, 250)) do
        fig = Figure();
        ax = Axis(fig[1, 1])
        graphplot!(ax, g; nlabels = labels)
        hidedecorations!(ax); hidespines!(ax)
        return fig
    end
end

# Some example data.
rng = MersenneTwister(1234)
n = 500
a1 = randn(rng, n)
a2 = randn(rng, n) .+ a1
a3 = randn(rng, n) .+ a2
a4 = randn(rng, n) .+ a2
a5 = randn(rng, n) .+ a3 .+ a4
x = [a1, a2, a3, a4, a5]

# Infer a directed graph using correlation-based independence tests
pairwise_test = SurrogateTest(PearsonCorrelation())
conditional_test = SurrogateTest(PartialCorrelation()) 
alg = PCRobust(pairwise_test, conditional_test; α = 0.5)
dg = infer_graph(alg, x)
```

Let's plot the resulting graph:

```@example causalgraph_corr
plotgraph(dg; labels = ["n$i" for i = 1:5])
```
