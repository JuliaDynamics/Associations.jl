using Test
using Graphs: SimpleDiGraph
using StableRNGs
using
rng = StableRNG(123)

function plotgraph(g)
    colors = [:black for i in 1:nv(g)]
    colors[1] = :red
    colors[2] = :green
    f, ax, p = graphplot(g,
        nlabels = repr.(1:nv(g)),
        nlabels_color = colors,
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

# Compare with CausalInference.jl
n = 10000
x = randn(rng, n)
y = x + 0.2*randn(rng, n)
z = x + 0.2*randn(rng, n)
w = y + z + 0.2*randn(rng, n)
q = w + 0.2*randn(rng, n)
r = w + 0.2*randn(rng, n)
X = [x, y, z, w, q, r]
df = (x=x, y=y, z=z, w=w, q=q, r=r)

alg = PC(CorrTest(), CorrTest(), α = 0.01)
sg, dg = infer_graph(alg, X; verbose = true)

# Test with different common independence tests (not exhaustive, because it takes too
# long for the tests to run then).
x, y, z = rand(rng, 70), rand(70), rand(70)
alg = PC(utest, ctest)

g = infer_graph(alg, x)
@test g isa SimpleDiGraph

# Example from CausalInference.jl docs
x = randn(rng, n)
v = x + randn(rng, n)*0.25
w = x + randn(rng, n)*0.25
z = v + w + randn(rng, n)*0.25
s = z + randn(rng, n)*0.25
X = [x, v, w, z, s]
df = (x=x, v=v, w=w, z=z, s=s)
