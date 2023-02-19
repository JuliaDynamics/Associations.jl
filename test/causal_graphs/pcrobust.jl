using Test
using Graphs: SimpleDiGraph
using Random; rng = Random.default_rng()

x = randn(rng, 100), randn(rng, 100), randn(rng, 100)
utest = SurrogateTest(MIShannon(), KSG1(); nshuffles = 4, rng)
ctest = SurrogateTest(PartialCorrelation(); nshuffles = 4, rng)
alg = PCRobust(utest, ctest)

g = infer_graph(alg, x)
@test g isa SimpleDiGraph
