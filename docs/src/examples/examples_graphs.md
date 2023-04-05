# Causal graphs

## [Optimal causation entropy](@id oce_example)

Here, we use the [`OCE`](@ref) algorithm to infer a time series graph. We use a
[`SurrogateTest`](@ref) for the initial step, and a [`LocalPermutationTest`](@ref)
for the conditional steps.

```@example causalgraph_oce
using CausalityTools
using Graphs
using Random
rng = MersenneTwister(1234)

# An example system where `X → Y → Z → W`.
sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(trajectory(sys, 400, Ttr = 10000))

# Independence tests for unconditional and conditional stages.
utest = SurrogateTest(MIShannon(), KSG2(k = 3, w = 1); rng, nshuffles = 150)
ctest = LocalPermutationTest(CMIShannon(), MesnerShalizi(k = 3, w = 1); rng, nshuffles = 150)

# Infer graph
alg = OCE(; utest, ctest, α = 0.05, τmax = 1)
parents = infer_graph(alg, [x, y, z, w])

# Convert to graph and inspect edges
g = SimpleDiGraph(parents)
collect(edges(g))
```

The algorithm nicely recovers the true causal directions.
