# Causal graphs

## [Optimal causation entropy](@id oce_example)

Here, we use the [`OCE`](@ref) algorithm to infer a time series graph. We use a
[`SurrogateTest`](@ref) for the initial step, and a [`LocalPermutationTest`](@ref)
for the conditional steps.

```@docs
using CausalityTools
using Random
rng = MersenneTwister(1234)

# An example system where `X → Y → Z → W`.
sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(trajectory(sys, 400, Ttr = 10000))

# Independence tests for unconditional and conditional stages.
utest = SurrogateTest(MIShannon(), KSG2(k = 5))
ctest = LocalPermutationTest(CMIShannon(), FPVP(k = 5))

# Infer graph
alg = OCE(; utest, ctest, α = 0.05, τmax = 3)
infer_graph(alg, [x, y, z, w])
```

The algorithm nicely recovers the true causal directions.

## [The PB-robust algorithm](@id pc_robust_example)

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
    with_theme(theme_minimal(), resolution = (400, 350)) do
        fig = Figure();
        ax = Axis(fig[1, 1])
        graphplot!(ax, g; nlabels = labels)
        hidedecorations!(ax); hidespines!(ax)
        return fig
    end
end

# Some example data.
rng = MersenneTwister(1234)

# The true graph is X → Y → Z → W
sys = system(Logistic4Chain(; rng))
X, Y, Z, W = columns(trajectory(sys, 1000, Ttr = 10000))
data = [X, Y, Z, W]

# Infer a directed graph using correlation-based independence tests
pairwise_test = SurrogateTest(TEShannon(), Zhu1(k = 10))
conditional_test = SurrogateTest(TEShannon(), Zhu1(k = 10)) 
alg = PCRobust(pairwise_test, conditional_test; α = 0.05)
g = infer_graph(alg, data)
```

Let's plot the resulting graph:

```@example causalgraph_corr
plotgraph(g; labels = ["a$i" for i = 1:5])
```
