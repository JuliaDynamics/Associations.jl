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
``x_2 \\to x_3`` and ``x_1 \\to x_3``. Modified from Gourévitch et al.
(2006)[Gourévitch2006].

## Equations of motion

```math
\\begin{aligned}
x_1(t+1) &= a_1 x_1 (1-x_1(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{1}(t) \\\\
x_2(t+1) &= a_1 x_2 (1-x_2(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{2}(t) + b x_1 x_2 \\\\
x_3(t+1) &= a_3 x_3 (1-x_3(t))^2  e^{-x_3(t)^2} + 0.4 \\xi_{3}(t) + c x_{2}(t) + d x_{1}(t)^2.
\\end{aligned}
```

[Gourévitch2006]:
    Gourévitch, B., Le Bouquin-Jeannès, R., & Faucon, G. (2006). Linear and nonlinear
    causality between signals: methods, examples and neurophysiological
    applications. Biological Cybernetics, 95(4), 349–369.
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
npts = 200
X, t = trajectory(sys, npts)


# Independence tests for unconditional and conditional stages.
uest = KSG2(MIShannon(); k=3, w=1)
utest = SurrogateAssociationTest(uest; rng, nshuffles=19)
cest = MesnerShalizi(CMIShannon(); k=3, w=1)
ctest = LocalPermutationTest(cest; rng, nshuffles=19)

# Infer graph
alg = PCMCI(; utest, ctest, α=0.05, τmax=1)
parents = infer_graph(alg, X)

# Convert to graph and inspect edges
g = SimpleDiGraph(parents)



using Graphs, GLMakie, GraphMakie

# Build edge labels aligned with `collect(edges(g))` (the order GraphMakie draws edges in),
# so each label lands on the correct edge. For every edge `src → dst`, we look up the
# `PCMCIParent`(s) in `res.parents[dst]` that are driven by `src`, and show their lag `τ`
# and test-statistic magnitude `|I|`. When several lags share the same `(src, dst)` pair
# (possible for `τmax > 1`, since a `SimpleDiGraph` collapses them into one edge), the
# labels are combined onto separate lines.
function pcmci_edge_labels(res::PCMCIResult, g)
    labels = String[]
    for e in edges(g)
        s, d = src(e), dst(e)
        ps = [p for p in res.parents[d] if p.i == s]
        if isempty(ps)
            push!(labels, "") # e.g. reverse half of a one-sided contemporaneous link
        else
            push!(labels,
                join(["τ=$(p.τ), |I|=$(round(abs(p.test_statistic), digits=2))" for p in ps], "\n"))
        end
    end
    return labels
end

function plotgraph(g; nlabels=repr.(1:nv(g)), elabels=nothing)
    f, ax, p = graphplot(g,
        ilabels=nlabels,
        ilabels_color=[:white for i in 1:nv(g)],
        node_color=:blue,
        node_size=80,
        arrow_size=15,
        elabels=elabels,
        elabels_color=:black,
        elabels_fontsize=11,
        curve_distance=0.2,       # separate the two arcs of a bidirectional (undirected) link
        curve_distance_usage=true,
        figure_padding=10
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
parents
elabels = pcmci_edge_labels(parents, g)
plotgraph(g, nlabels=["x1", "x2", "x3"], elabels=elabels)


# ---------------------------------------------------------------------------
# Five variables with longer lags
# ---------------------------------------------------------------------------
# A 5-variable linear system where the couplings act at different (and longer)
# lags. The ground-truth directed links are:
#
#   x1(t-2) → x2(t)
#   x2(t-3) → x3(t)
#   x1(t-4) → x4(t)
#   x3(t-1) → x5(t)
#   x4(t-5) → x5(t)
#
# In addition, every variable has an autoregressive dependence on its own
# previous value. Because the largest interaction lag is 5, we need `τmax ≥ 5`
# for PCMCI to be able to recover all of these links.

using StateSpaceSets: StateSpaceSet

function linear5(npts::Int; rng=Random.default_rng(),
    a=0.5,      # autoregressive strength
    c=0.5,      # coupling strength
    σ=0.4)      # noise strength
    # `pad` extra transient points are generated (and later discarded) so that
    # the longest lag (5) always has valid history.
    pad = 100
    N = npts + pad
    x1 = zeros(N);
    x2 = zeros(N);
    x3 = zeros(N);
    x4 = zeros(N);
    x5 = zeros(N)
    for t in 6:N
        x1[t] = a*x1[t-1] + σ*randn(rng)
        x2[t] = a*x2[t-1] + c*x1[t-2] + σ*randn(rng)
        x3[t] = a*x3[t-1] + c*x2[t-3] + σ*randn(rng)
        x4[t] = a*x4[t-1] + c*x1[t-4] + σ*randn(rng)
        x5[t] = a*x5[t-1] + c*x3[t-1] + c*x4[t-5] + σ*randn(rng)
    end
    return StateSpaceSet(x1[(pad+1):end], x2[(pad+1):end], x3[(pad+1):end],
        x4[(pad+1):end], x5[(pad+1):end])
end

rng5 = Xoshiro(1234)
X5 = linear5(200; rng=rng5)

# Linear association measures are appropriate here, and cheaper than the
# nearest-neighbour estimators used above.
uest = KSG2(MIShannon(); k=3, w=1)
utest = SurrogateAssociationTest(uest; rng, nshuffles=19)
cest = MesnerShalizi(CMIShannon(); k=3, w=1)
ctest = LocalPermutationTest(cest; rng, nshuffles=19)

# `τmax = 5` so the longest coupling (x4 → x5 at lag 5) can be detected.
alg5 = PCMCI(; utest=utest, ctest=ctest, α=0.05, τmax=5)
parents5 = infer_graph(alg5, X5, verbose=true)

g5 = SimpleDiGraph(parents5)
elabels5 = pcmci_edge_labels(parents5, g5)
plotgraph(g5, nlabels=["x1", "x2", "x3", "x4", "x5"], elabels=elabels5)