using Test
using Random
rng = Xoshiro(1234)

# We can use surrogate tests and p-values to further verify the correctness of the 
# algorithm.
test = SurrogateAssociationTest(ChatterjeeCorrelation(), nshuffles = 19)

# We expect that we *cannot* reject the null hypothesis for independent variables
x = rand(rng, 1:10, 100)
y = rand(rng, 1:10, 100)
α = 0.05
@test pvalue(independence(test, x, y)) > α 

# We expect that we *can* reject the null hypothesis for for dependent variables.
w = rand(rng, 1:10, 100)
z = rand(rng, 1:10, 100) .* sin.(w) .+ cos.(w)
@test pvalue(independence(test, z, w)) < α 
