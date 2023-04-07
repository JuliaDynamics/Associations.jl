using Test
using Graphs: SimpleDiGraph
using StableRNGs
rng = StableRNG(123)

# Compare with CausalInference.jl
n = 5000
x = randn(rng, n)
y = x + 0.25*randn(rng, n)
z = x + 0.25*randn(rng, n)
w = y + z + 0.25*randn(rng, n)
q = w + 0.25*randn(rng, n)
r = w + 0.25*randn(rng, n)
X = [x, y, z, w, q, r]
alg = PC(CorrTest(), CorrTest())

sg, dg = infer_graph(alg, X; verbose = true)

# Test with different common independence tests (not exhaustive, because it takes too
# long for the tests to run then).
x, y, z = rand(rng, 70), rand(70), rand(70)
alg = PC(utest, ctest)

g = infer_graph(alg, x)
@test g isa SimpleDiGraph

function plotgraph(g)
    colors = [:black for i in 1:nv(g)]
    colors[1] = :red
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
