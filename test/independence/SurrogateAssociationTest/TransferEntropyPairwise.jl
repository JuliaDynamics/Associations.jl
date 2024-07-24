using Random
rng = Xoshiro(1234)

sys = system(Logistic2Unidir(; c_xy = 0.5))
x, y = columns(first(trajectory(sys, 300, Ttr = 10000)))

# Creation
@test SurrogateAssociationTest(est) isa SurrogateAssociationTest
est = CMIDecomposition(TEShannon(), FPVP())
# ArgumentError thrown if an estimator isn't provided.
@test_throws ArgumentError SurrogateAssociationTest(TEShannon())


α = 0.0287 # Arbitrary significance level 1 - α = 0.9713
test = SurrogateAssociationTest(est; rng, nshuffles = 50)

# The ground truth is X → Y, so we should be able to reject the null
# when testing transferentropy(x → y)
@test pvalue(independence(test, x, y)) < α

# The ground truth is X → Y, so we shouldn't be able to reject the null
# when testing transferentropy(y → x)
@test pvalue(independence(test, y, x)) > α

x, y = columns(first(trajectory(sys, 100, Ttr = 1000)))
est = Lindner()
test = SurrogateAssociationTest(est, nshuffles = 19)
@test independence(test, x, y) isa SurrogateAssociationTestResult

est = Zhu1()
test = SurrogateAssociationTest(est; nshuffles = 19)
@test independence(test, x, y) isa SurrogateAssociationTestResult
