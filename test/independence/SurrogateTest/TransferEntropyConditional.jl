# A chain of coupled logistic maps. Set the coupling from first to second
# variable high, so that transferentropy(x → z) becomes significant.
# This should vanish when doing transferentropy(x → z | y)
using Random
rng = Random.MersenneTwister(1234)
sys = system(Logistic4Chain(; xi = [0.1, 0.2, 0.3, 0.4], rng));
n = 500
x, y, z, w = columns(first(trajectory(sys, n, Ttr = 10000)));

α = 0.04 # Arbitrary significance level 1 - α = 0.96

# The ground truth is X → Y → Z.
test = SurrogateTest(TEShannon(), FPVP(); rng)


@test pvalue(independence(test, x, z)) < α # This has been manually tested to occur with c₁₂ = 0.8

# We should be able to reject the null when testing transferentropy(x → y | z)
@test pvalue(independence(test, x, z, y)) > α

@test independence(SurrogateTest(TEShannon(), Zhu1()), x, y, z) isa SurrogateTestResult
@test independence(SurrogateTest(TEShannon(), Lindner()), x, y, z) isa SurrogateTestResult
