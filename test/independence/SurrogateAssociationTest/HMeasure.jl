# Analytical tests (in the limit of a lot of "many" samples)
# ------------------------------------------------------------
using Random
rng = MersenneTwister(1234)
n = 100
x, y = rand(rng, n), rand(rng, n)
z = x .+ y
test = SurrogateAssociationTest(HMeasure(); rng)
α = 0.04 # Some arbitrary significance level.

# We shouldn't be able to reject the null when the variables are independent
@test pvalue(independence(test, x, y)) > 0.04

# We should reject the null when variables are dependent.
@test pvalue(independence(test, x, z)) < 0.04
