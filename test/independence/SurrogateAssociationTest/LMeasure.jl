# Analytical tests (in the limit of a lot of samples)
# ------------------------------------------------------------
using Random
rng = StableRNG(1234)
n = 100
x, y = rand(rng, n), rand(rng, n)
z = x .+ y
test = SurrogateAssociationTest(LMeasure(); rng)
α = 0.05 # Some arbitrary significance level.

# We shouldn't be able to reject the null when the variables are independent
@test pvalue(independence(test, x, y)) > α

# We should reject the null when variables are dependent.
@test pvalue(independence(test, x, z)) < α
