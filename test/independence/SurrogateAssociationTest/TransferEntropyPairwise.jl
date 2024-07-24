using Test
using Random
using TimeseriesSurrogates

rng = Xoshiro(1234)

sys = system(Logistic2Unidir(; c_xy = 0.5, rng))
x, y = columns(first(trajectory(sys, 400, Ttr = 10000)))

# Creation
est = MIDecomposition(TEShannon(), KSG1())
@test SurrogateAssociationTest(est) isa SurrogateAssociationTest
est = CMIDecomposition(TEShannon(), FPVP())
# ArgumentError thrown if an estimator isn't provided.
@test_throws ArgumentError SurrogateAssociationTest(TEShannon())


α = 0.04 # Arbitrary significance level 1 - α = 0.96
test = SurrogateAssociationTest(est; rng, nshuffles = 19)

# The ground truth is X → Y, so we should be able to reject the null
# when testing transferentropy(x → y)
@test pvalue(independence(test, x, y)) < α

# The ground truth is X → Y, so we shouldn't be able to reject the null
# when testing transferentropy(y → x)
@test pvalue(independence(test, y, x)) > α

# MIDecomposition estimator 
est = MIDecomposition(TEShannon(), KSG1())
test = SurrogateAssociationTest(est, nshuffles = 2)
@test independence(test, x, y) isa SurrogateAssociationTestResult

# Dedicated estimators
x, y = columns(first(trajectory(sys, 200, Ttr = 1000)))
est = Lindner()
test = SurrogateAssociationTest(est, nshuffles = 2)
@test independence(test, x, y) isa SurrogateAssociationTestResult

est = Zhu1()
test = SurrogateAssociationTest(est; nshuffles = 2)
@test independence(test, x, y) isa SurrogateAssociationTestResult

# `EntropyDecomposition`
est = EntropyDecomposition(TEShannon(), Kraskov())
test = SurrogateAssociationTest(est; nshuffles = 2);
@test independence(test, x, y) isa SurrogateAssociationTestResult

# Can't use single-variable surrogate methods when a dimension is higher than 1
est = EntropyDecomposition(TEShannon(embedding = EmbeddingTE(dS = 2)), Kraskov())
test = SurrogateAssociationTest(est; nshuffles = 2, surrogate = AAFT());
@test_throws ArgumentError independence(test, x, y)

# Optimising parameters using traditional methods
est = EntropyDecomposition(TEShannon(embedding = OptimiseTraditional()), Kraskov())
test = SurrogateAssociationTest(est; nshuffles = 2);
@test independence(test, x, y) isa SurrogateAssociationTestResult
