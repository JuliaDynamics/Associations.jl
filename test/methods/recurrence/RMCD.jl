using StableRNGs
using StateSpaceSets
using Distances: Chebyshev

rng = StableRNG(1234)
n = 100
x = rand(rng, n)
y = rand(rng, n)
z = rand(rng, n)
X = rand(rng, n, 2) |> StateSpaceSet
Y = rand(rng, n, 3) |> StateSpaceSet
Z = rand(rng, n, 2) |> StateSpaceSet

@test_throws UndefKeywordError RMCD()

@test association(RMCD(; r = 0.5), x, y) >= 0
@test association(RMCD(; r = 0.5), x, Y) >= 0
@test association(RMCD(; r = 0.5), X, Y) >= 0
@test association(RMCD(; r = 0.5, metric = Chebyshev()), x, y) >= 0
@test association(RMCD(; r = 0.5, metric = Chebyshev()), X, Y) >= 0

@test association(RMCD(; r = 0.5), x, y, z) >= 0
@test association(RMCD(; r = 0.5), x, Y, z) >= 0
@test association(RMCD(; r = 0.5), X, Y, Z) >= 0
@test association(RMCD(; r = 0.5, metric = Chebyshev()), x, y, z) >= 0
@test association(RMCD(; r = 0.1, metric = Chebyshev()), X, Y, z) >= 0
@test association(RMCD(; r = 0.5), x, y, x) == 0
@test association(RMCD(; r = 0.5), x, y, y) == 0

# We should not be able to reject null for independent variables
test = SurrogateAssociationTest(RMCD(r = 0.5); rng, nshuffles = 50)

@test pvalue(independence(test, x, y)) >= α
@test pvalue(independence(test, X, Y)) >= α
@test pvalue(independence(test, x, Y)) >= α

# Test on a dynamical system.
sys = system(Logistic4Chain(; xi = rand(rng, 4), rng))
x, y, z, w = columns(first(trajectory(sys, 200, Ttr = 10000)));
test = LocalPermutationTest(RMCD(r = 0.5); rng)

# X and Z are independent given Y, so we shouldn't be able to reject the null (p > α).
pval = pvalue(independence(test, x, z, y))
α = 0.05
@test pval > α
