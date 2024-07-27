using Random: MersenneTwister
rng = MersenneTwister(12346)
x, y = randn(rng, 1000), randn(rng, 1000)
m = JointDistanceDistribution(D = 3, B = 5)
test =  JointDistanceDistributionTest(m)
@test test isa JointDistanceDistributionTest
@test independence(test, x, y) isa Associations.JDDTestResult

# Don't reject null at significance level (1 - α) when there is no coupling.
α = 0.05
@test pvalue(independence(test, x, y)) > α 

# Reject null at significance level (1 - α) when there is coupling
z = y .+ x
@test pvalue(independence(test, z, x)) < α

# Test summary outouts
out_str = repr(independence(test, x, y))
occursin("Independence cannot be rejected", out_str)