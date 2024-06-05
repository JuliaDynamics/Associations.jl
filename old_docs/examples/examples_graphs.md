# Causal graphs

Before introducing the causal graph examples, let's create a function that can plot
directed graphs that we'll use below.

```@example graph_examples
using Graphs, CairoMakie, GraphMakie

function plotgraph(g)
    f, ax, p = graphplot(g,
        nlabels = repr.(1:nv(g)),
        nlabels_color = [:red for i in 1:nv(g)],
    )
    offsets = 0.05 * (p[:node_pos][] .- p[:node_pos][][1])
    offsets[1] = Point2f(0, 0.2)
    p.nlabels_offset[] = offsets
    autolimits!(ax)
    hidedecorations!(ax)
    hidespines!(ax)
    ax.aspect = DataAspect()
    return f
end
```

## [Optimal causation entropy](@id oce_example)

Here, we use the [`OCE`](@ref) algorithm to infer a time series graph. We use a
[`SurrogateTest`](@ref) for the initial step, and a [`LocalPermutationTest`](@ref)
for the conditional steps.

```@example graph_examples
using CausalityTools
using StableRNGs
rng = StableRNG(123)

# An example system where `X → Y → Z → W`.
sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(first(trajectory(sys, 400, Ttr = 10000)))

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

The algorithm nicely recovers the true causal directions. We can also plot the graph using
the function we made above.

```@example graph_examples
plotgraph(g)
```

## [PC-algorithm](@id pc_examples)

### [Correlation-based tests](@id pc_examples_corr)

Here, we demonstrate the use of the [`PC`](@ref)-algorithm with the correlation-based
[`CorrTest`](@ref) both for the pairwise (i.e. using [`PearsonCorrelation`](@ref))
and conditional (i.e. using [`PartialCorrelation`](@ref)) case.

We'll reproduce the first example from
[CausalInference.jl](https://github.com/mschauer/CausalInference.jl), where they
also use a parametric correlation test to infer the skeleton graph for some 
normally distributed data.

```@example graph_examples
using CausalityTools
using StableRNGs
rng = StableRNG(123)
n = 500
v = randn(rng, n)
x = v + randn(rng, n)*0.25
w = x + randn(rng, n)*0.25
z = v + w + randn(rng, n)*0.25
s = z + randn(rng, n)*0.25
X = [x, v, w, z, s]

# Infer a completed partially directed acyclic graph (CPDAG)
alg = PC(CorrTest(), CorrTest(); α = 0.05)
est_cpdag_parametric = infer_graph(alg, X; verbose = false)

# Plot the graph
plotgraph(est_cpdag_parametric)
```

### [Nonparametric tests](@id pc_examples_nonparametric)

The main difference between the PC algorithm implementation here and in
CausalInference.jl is that our implementation automatically works with any compatible
and [`IndependenceTest`](@ref), and thus any combination of (nondirectional)
[`AssociationMeasure`](@ref) and estimator.

Here, we replicate the example above, but using a nonparametric [`SurrogateTest`](@ref)
with the Shannon mutual information [`MIShannon`](@ref) measure and the
[`GaoOhViswanath`](@ref) estimator for the pairwise independence tests, and a
[`LocalPermutationTest`](@ref) with conditional mutual information [`CMIShannon`](@ref)
and the [`MesnerShalizi`](@ref).

```@example graph_examples
rng = StableRNG(123)

# Use fewer observations, because MI/CMI takes longer to estimate
n = 400
v = randn(rng, n)
x = v + randn(rng, n)*0.25
w = x + randn(rng, n)*0.25
z = v + w + randn(rng, n)*0.25
s = z + randn(rng, n)*0.25
X = [x, v, w, z, s]

pairwise_test = SurrogateTest(MIShannon(), GaoOhViswanath(k = 10))
cond_test = LocalPermutationTest(CMIShannon(), MesnerShalizi(k = 10))
alg = PC(pairwise_test, cond_test; α = 0.05)
est_cpdag_nonparametric = infer_graph(alg, X; verbose = false)
plotgraph(est_cpdag_nonparametric)
```

We get the same graph as with the parametric estimator. However, for general non-gaussian
data, the correlation-based tests (which assumes normally distributed data)
will *not* give the same results as other independence tests.
