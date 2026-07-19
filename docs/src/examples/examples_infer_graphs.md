# [Causal graphs](@id examples_network_inference)

Before introducing the causal graph examples, let's create a function that can plot
directed graphs that we'll use below.

```@example graph_examples
using Graphs, CairoMakie, GraphMakie

function plotgraph(g; nlabels = repr.(1:nv(g)))
    f, ax, p = graphplot(g,
        ilabels = nlabels,
        ilabels_color = [:white for i in 1:nv(g)],
        node_color = :blue,
        node_size = 80,
        arrow_size = 15,
        figure_padding = 10
    )
    offsets = 0.02 * (p[:node_pos][] .- p[:node_pos][][1])
    offsets[1] = Point2f(0, 0.2)
    p.nlabels_offset[] = offsets
    autolimits!(ax)
    hidedecorations!(ax)
    hidespines!(ax)
    ax.aspect = DataAspect()
    return f
end
```

We'll also implement a set of chained logistic maps with unidirectional coupling.

```@example graph_examples
using DynamicalSystemsBase
Base.@kwdef struct Logistic4Chain{V, RX, RY, RZ, RW, C1, C2, C3, Σ1, Σ2, Σ3, RNG}
    xi::V = [0.1, 0.2, 0.3, 0.4]
    rx::RX = 3.9
    ry::RY = 3.6
    rz::RZ = 3.6
    rw::RW = 3.8
    c_xy::C1 = 0.4
    c_yz::C2 = 0.4
    c_zw::C3 = 0.35
    σ_xy::Σ1 = 0.05
    σ_yz::Σ2 = 0.05
    σ_zw::Σ3 = 0.05
    rng::RNG = Random.default_rng()
end

function eom_logistic4_chain(u, p::Logistic4Chain, t)
    (; xi, rx, ry, rz, rw, c_xy, c_yz, c_zw, σ_xy, σ_yz, σ_zw, rng) = p
    x, y, z, w = u
    f_xy = (y +  c_xy*(x + σ_xy * rand(rng)) ) / (1 + c_xy*(1+σ_xy))
    f_yz = (z +  c_yz*(y + σ_yz * rand(rng)) ) / (1 + c_yz*(1+σ_yz))
    f_zw = (w +  c_zw*(z + σ_zw * rand(rng)) ) / (1 + c_zw*(1+σ_zw))
    dx = rx * x * (1 - x)
    dy = ry * (f_xy) * (1 - f_xy)
    dz = rz * (f_yz) * (1 - f_yz)
    dw = rw * (f_zw) * (1 - f_zw)
    return SVector{4}(dx, dy, dz, dw)
end


function system(definition::Logistic4Chain)
    return DiscreteDynamicalSystem(eom_logistic4_chain, definition.xi, definition)
end
```


## [Optimal causation entropy](@id oce_example)

Here, we use the [`OCE`](@ref) algorithm to infer a time series graph. We use a
[`SurrogateAssociationTest`](@ref) for the initial step, and a [`LocalPermutationTest`](@ref)
for the conditional steps.

```@example graph_examples
using Associations
using StableRNGs
rng = StableRNG(123)

# An example system where `X → Y → Z → W`.
sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(first(trajectory(sys, 300, Ttr = 10000)))

# Independence tests for unconditional and conditional stages.
uest = KSG2(MIShannon(); k = 3, w = 1)
utest = SurrogateAssociationTest(uest; rng, nshuffles = 19)
cest =  MesnerShalizi(CMIShannon(); k = 3, w = 1)
ctest = LocalPermutationTest(cest; rng, nshuffles = 19)

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
plotgraph(g; nlabels = ["x", "y", "z", "w"])
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
using Associations
using StableRNGs
rng = StableRNG(123)
n = 300
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
plotgraph(est_cpdag_parametric; nlabels = ["x", "v", "w", "z", "s"])
```

### [Nonparametric tests](@id pc_examples_nonparametric)

The main difference between the PC algorithm implementation here and in
CausalInference.jl is that our implementation automatically works with any compatible
and [`IndependenceTest`](@ref), and thus any combination of (nondirectional)
[`AssociationMeasure`](@ref) and estimator.

Here, we replicate the example above, but using a nonparametric [`SurrogateAssociationTest`](@ref)
with the Shannon mutual information [`MIShannon`](@ref) measure and the
[`GaoOhViswanath`](@ref) estimator for the pairwise independence tests, and a
[`LocalPermutationTest`](@ref) with conditional mutual information [`CMIShannon`](@ref)
and the [`MesnerShalizi`](@ref).

```@example graph_examples
using Associations
using StableRNGs
rng = StableRNG(123)

# Use fewer observations, because MI/CMI takes longer to estimate
n = 300
v = randn(rng, n)
x = v + randn(rng, n)*0.25
w = x + randn(rng, n)*0.25
z = v + w + randn(rng, n)*0.25
s = z + randn(rng, n)*0.25
X = [x, v, w, z, s]

est_pairwise = JointProbabilities(MIShannon(), CodifyVariables(ValueBinning(3)))
est_cond = MesnerShalizi(CMIShannon(); k = 5)
pairwise_test = SurrogateAssociationTest(est_pairwise; rng, nshuffles = 50)
cond_test = LocalPermutationTest(est_cond; rng, nshuffles = 50)
alg = PC(pairwise_test, cond_test; α = 0.05)
est_cpdag_nonparametric = infer_graph(alg, X; verbose = false)
plotgraph(est_cpdag_nonparametric; nlabels = ["x", "v", "w", "z", "s"])
```

We get the same basic structure of the graph, but which directional associations 
are correctly ruled out varies. In general, using different types of 
association measures with different independence tests, applied to general 
non-gaussian data, will not give the same results as the correlation-based tests.

## [PCMCI-algorithm](@id pcmci_examples)

```@example graph_examples_pcmci
using Associations
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem, trajectory
using Distributions: Normal
using Random

export Nonlinear3
"""
    Nonlinear3 <: DiscreteDefinition
    Nonlinear3(; xi = rand(3),
        σ₁ = 1.0, σ₂ = 1.0, σ₃ = 1.0,
        a₁ = 3.4, a₂ = 3.4, a₃ = 3.4,
        b₁ = 0.4, b₂ = 0.4, b₃ = 0.4,
        c₁₂ = 0.5, c₂₃ = 0.3, c₁₃ = 0.5,
        rng = Random.default_rng())

A 3d nonlinear system with nonlinear couplings ``x_1 \\to x_2``,
``x_2 \\to x_3`` and ``x_1 \\to x_3``. Modified from Gourévitch, B., Le Bouquin-Jeannès, R., & Faucon, G. (2006). Linear and nonlinear causality between signals: methods, examples and neurophysiological applications. Biological Cybernetics, 95(4), 349–369.
"""
Base.@kwdef struct Nonlinear3{V,Σx,Σy,Σz,AX,AY,AZ,BX,BY,BZ,C1,C2,C3,RNG}
    xi::V = [0.1, 0.2, 0.3]
    σx::Σx = Normal(0, 1.0)
    σy::Σy = Normal(0, 1.0)
    σz::Σz = Normal(0, 1.0)
    ax::AX = 3.4
    ay::AY = 3.4
    az::AZ = 3.4
    bx::BX = 0.4
    by::BY = 0.4
    bz::BZ = 0.4
    c_xy::C1 = 0.5
    c_xz::C2 = 0.3
    c_yz::C3 = 0.5
    rng::RNG = Random.default_rng()
end

function system(definition::Nonlinear3)
    return DiscreteDynamicalSystem(eom_nonlinear3, definition.xi, definition)
end

function eom_nonlinear3(u, p, n)
    x, y, z = u
    (; xi, σx, σy, σz, ax, ay, az, bx, by, bz, c_xy, c_xz, c_yz, rng) = p
    ξ₁ = rand(rng, σx)
    ξ₂ = rand(rng, σy)
    ξ₃ = rand(rng, σz)
    dx = ax * x * (1-x)^2 * exp(-x^2) + bx*ξ₁
    dy = ay * y * (1-y)^2 * exp(-y^2) + by*ξ₂ + c_xy*x*y
    dz = az * z * (1-z)^2 * exp(-z^2) + bz*ξ₃ + c_yz*y + c_xz*x^2
    return SVector{3}(dx, dy, dz)
end

rng = Xoshiro(1234)
sys = system(Nonlinear3(; rng))
npts = 100
X, t = trajectory(sys, npts)

# Independence tests for unconditional and conditional stages.
uest = KSG2(MIShannon(); k=3, w=1)
utest = SurrogateAssociationTest(uest; rng, nshuffles=19)
cest = MesnerShalizi(CMIShannon(); k=3, w=1)
ctest = LocalPermutationTest(cest; rng, nshuffles=19)

# Infer graph
alg = PCMCI(; utest, ctest, α=0.05, τmax=1)
est = infer_graph(alg, X; verbose=true)

plotgraph(est; nlabels = ["x", "y", "z"])
```
