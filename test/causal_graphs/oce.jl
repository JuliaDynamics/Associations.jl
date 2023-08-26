using CausalityTools
using CausalityTools: OCESelectedParents
using Test
using StableRNGs
using Graphs.SimpleGraphs: SimpleEdge

rng = StableRNG(123)
sys = system(Logistic4Chain(; rng))
X = columns(first(trajectory(sys, 60, Ttr = 10000)))
utest = SurrogateTest(MIShannon(), KSG1(k = 10, w = 1); rng, nshuffles = 30)
ctest = LocalPermutationTest(CMIShannon(), MesnerShalizi(k = 10, w = 1); rng, nshuffles = 30)
alg = OCE(; utest, ctest, τmax = 2)
parents = infer_graph(alg, X; verbose = true)
@test parents isa Vector{<:OCESelectedParents}
@test SimpleDiGraph(parents) isa SimpleDiGraph

rng = StableRNG(123)
sys = system(Logistic2Bidir(; rng))
X = columns(first(trajectory(sys, 300, Ttr = 10000)))
utest = SurrogateTest(MIShannon(), KSG1(k = 10, w = 1); rng, nshuffles = 100)
ctest = LocalPermutationTest(CMIShannon(), MesnerShalizi(k = 10, w = 1); rng, nshuffles = 100)
parents = infer_graph(OCE(; utest, ctest, τmax = 1), X; verbose = true)
@test parents isa Vector{<:OCESelectedParents}
g = SimpleDiGraph(parents)
@test g isa SimpleDiGraph

# "Analytical" test: check that we at least identify one true positive. There may
# be several false positives, but there's no way of telling a priori how many.
function at_least_one_true_positive(true_edges, estimated_graph)
    estimated_edges = edges(estimated_graph)
    at_least_one_tp = false
    for e in true_edges
        if e in estimated_edges
            at_least_one_tp = true
        end
    end
    return at_least_one_tp
end

@test at_least_one_true_positive([SimpleEdge(1, 2), SimpleEdge(2, 1)], g)
