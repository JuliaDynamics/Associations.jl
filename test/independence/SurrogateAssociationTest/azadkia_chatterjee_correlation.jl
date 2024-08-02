using Test
using Random
rng = Xoshiro(1234)

# We can use surrogate tests and p-values to further verify the correctness of the 
# algorithm.
test = SurrogateAssociationTest(AzadkiaChatterjeeCoefficient(), nshuffles = 19, rng = rng)
n = 200
# We expect that we *cannot* reject the null hypothesis for independent variables
x = rand(rng, n)
y = rand(rng, n)
α = 0.05
@test pvalue(independence(test, x, y)) > α 

# We expect that we *can* reject the null hypothesis for for dependent variables.
x = rand(rng, n)
y = rand(rng, n) .* x
@test pvalue(independence(test, x, y )) < α 

# We expect that we *cannot* reject the null hypothesis for two extremal variables 
# connected by an intermediate variable when conditioning on the intermediate variable.
x = rand(rng, n)
y = rand(rng, n) .+ x
z = rand(rng, n) .* y
@test pvalue(independence(test, x, z, y)) > α 
