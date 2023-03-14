using StableRNGs
using StateSpaceSets
using Distances: Chebyshev

rng = StableRNG(1234)
n = 200
x = rand(rng, n)
y = rand(rng, n)
z = rand(rng, n)
X = rand(rng, n, 2) |> StateSpaceSet
Y = rand(rng, n, 3) |> StateSpaceSet
Z = rand(rng, n, 2) |> StateSpaceSet

@test_throws UndefKeywordError RMCD()

@test rmcd(RMCD(; r = 0.5), x, y) >= 0
@test rmcd(RMCD(; r = 0.5), x, Y) >= 0
@test rmcd(RMCD(; r = 0.5), X, Y) >= 0
@test rmcd(RMCD(; r = 0.5, metric = Chebyshev()), x, y) >= 0
@test rmcd(RMCD(; r = 0.5, metric = Chebyshev()), X, Y) >= 0

@test rmcd(RMCD(; r = 0.5), x, y, z) >= 0
@test rmcd(RMCD(; r = 0.5), x, Y, z) >= 0
@test rmcd(RMCD(; r = 0.5), X, Y, Z) >= 0
@test rmcd(RMCD(; r = 0.5, metric = Chebyshev()), x, y, z) >= 0
@test rmcd(RMCD(; r = 0.1, metric = Chebyshev()), X, Y, z) >= 0
@test rmcd(RMCD(; r = 0.5), x, y, x) == 0
@test rmcd(RMCD(; r = 0.5), x, y, y) == 0

test = SurrogateTest(RMCD(r = 0.2); rng)
α = 0.05
x, y, z, w = columns(trajectory(system(Logistic4Chain(; xi = rand(4), rng)), 500, Ttr = 10000));
# We should not be able to reject null for independent variables
@test pvalue(independence(test, x, y)) >= α
@test pvalue(independence(test, X, Y)) >= α
@test pvalue(independence(test, x, Y)) >= α

# For dependent variables, we should be able to reject the null hypothesis of
# independence. This goes both ways.
rng = StableRNG(1234)
n = 500
x = rand(rng, n)
z = y .+ rand(rng, n)
@test pvalue(independence(test, x, z)) < α
rng = StableRNG(1234)
test = SurrogateTest(RMCD(r = 0.2); rng)
@test pvalue(independence(test, z, x)) < α
