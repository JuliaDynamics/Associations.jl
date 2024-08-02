using Test
using Random
rng = Xoshiro(1234)
n = 200
# We expect that we *cannot* reject the null hypothesis for two extremal variables 
# connected by an intermediate variable when conditioning on the intermediate variable.
test = LocalPermutationTest(AzadkiaChatterjeeCoefficient(), nshuffles = 19, rng = rng)
x = rand(rng, n)
y = rand(rng, n) .+ x
z = rand(rng, n) .* y
@test pvalue(independence(test, x, z, y)) > α 
