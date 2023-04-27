# Analytical tests (in the limit of a lot of samples)
# ------------------------------------------------------------
using Random
rng = MersenneTwister(1234)
x, y = rand(rng, 500), rand(rng, 500)
z = x .+ y
test = SurrogateTest(LMeasure(); rng)
α = 0.05 # Some arbitrary significance level.

# We shouldn't be able to reject the null when the variables are independent
@test pvalue(independence(test, x, y)) > α

# We should reject the null when variables are dependent.
@test pvalue(independence(test, x, z)) < α
