using Associations
using Associations: PCMCIResult
using DynamicalSystemsBase
using Test
using StableRNGs

# Test on a real simulated system where X1 -> X2 -> X3 -> X4
# with lag -1 for each link.
# (we don't expect this to be 100% successful in general, but this example works)
rng = StableRNG(1234);
sys = system(Logistic4Chain(; rng));
X = first(trajectory(sys, 200, Ttr = 10000));
utest = SurrogateAssociationTest(KSG1(k = 10), nshuffles = 19);
ctest = SurrogateAssociationTest(Rahimzamani(k = 10), nshuffles = 19);
alg = PCMCI(; utest, ctest, qmax = 1, pmax = 5, τmax = 1);

res = infer_graph(alg, X)

@test res isa PCMCIResult
@test isempty(res.parents[1])
@test length(res.parents[2]) == 1 # only one driver for X2
@test first(res.parents[2]).i == 1 && first(res.parents[2]).τ == -1 # the driver for X2(0) is X1(-1)
@test length(res.parents[3]) == 1 # only one driver for X3
@test first(res.parents[3]).i == 2 && first(res.parents[3]).τ == -1 # the driver for X2(0) is X2(-1)
@test length(res.parents[4]) == 1 # only one driver for X4
@test first(res.parents[4]).i == 3 && first(res.parents[4]).τ == -1 # the driver for X3(0) is X3(-1)

alg = PCMCI(; utest, ctest, qmax = 1, pmax = 5, τmax = 1, fdr_adjust = false);
res = infer_graph(alg, X)
@test res isa PCMCIResult