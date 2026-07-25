using Associations
using Associations: PCMCIResult
using DynamicalSystemsBase
using Test
using StableRNGs


# Test on a real simulated system where X1 -> X2 -> X3 -> X4
# with lag -1 for each link.
# (we don't expect this to be 100% successful in general, but this example works)
rng = StableRNG(12354);
sys = system(Logistic4Chain(; rng));
X = first(trajectory(sys, 500, Ttr=10000));
utest = SurrogateAssociationTest(ChatterjeeCorrelation(), nshuffles=19);
ctest = LocalPermutationTest(AzadkiaChatterjeeCoefficient(), nshuffles=19);

alg = PCMCI(;
    utest,
    ctest,
    qmax=2, pmax=3, τmax=2,
    use_abs_utest=true,
    use_abs_ctest=true);

res = infer_graph(alg, X)

# We only test that the algorithm runs. For this example, we don't expect 
# PCMCI to yield only the correct links, because it is a system of the type 
# the method is relatively poor at. 
@test res isa PCMCIResult
