# Analytical tests (in the limit of a lot of samples)
# ------------------------------------------------------------
using Random
rng = MersenneTwister(1234)
x, y = rand(rng, 300), rand(rng, 300)
z = x .+ y
test = SurrogateTest(SMeasure())
α = 0.04 # Some arbitrary significance level.

# We shouldn't be able to reject the null when the variables are independent
@test pvalue(independence(test, x, y)) > 0.04

# We should reject the null when variables are dependent.
@test pvalue(independence(test, x, z)) < 0.04
