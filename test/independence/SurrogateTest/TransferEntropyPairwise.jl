sys = logistic2_unidir(c_xy = 0.5)
x, y = columns(trajectory(sys, 3000, Ttr = 10000))

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