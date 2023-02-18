using Random
rng = Random.MersenneTwister(1234)
sys = system(Logistic2Unidir(; c_xy = 0.5, rng))
x, y = columns(trajectory(sys, 1000, Ttr = 10000))

# ArgumentError thrown if an estimator isn't provided.
@test_throws ArgumentError SurrogateTest(TEShannon())
@test SurrogateTest(TEShannon(), FPVP()) isa SurrogateTest

α = 0.0287 # Arbitrary significance level 1 - α = 0.9713
test = SurrogateTest(TEShannon(), FPVP())

# The ground truth is X → Y, so we should be able to reject the null
# when testing transferentropy(x → y)
@test pvalue(independence(test, x, y)) < α

# The ground truth is X → Y, so we shouldn't be able to reject the null
# when testing transferentropy(y → x)
@test pvalue(independence(test, y, x)) > α

x, y = columns(trajectory(sys, 100, Ttr = 1000))
@test independence(SurrogateTest(TEShannon(), Lindner()), x, y) isa SurrogateTestResult
@test independence(SurrogateTest(TEShannon(), Zhu1()), x, y) isa SurrogateTestResult
