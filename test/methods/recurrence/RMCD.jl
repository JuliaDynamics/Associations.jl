using StableRNGs
using StateSpaceSets
using Distances: Chebyshev

rng = StableRNG(1234)
n = 250
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

# We should not be able to reject null for independent variables
@test pvalue(independence(test, x, y)) >= α
@test pvalue(independence(test, X, Y)) >= α
@test pvalue(independence(test, x, Y)) >= α


# Test on a dynamical system.
sys = system(Logistic4Chain(; xi = rand(rng, 4), rng))
x, y, z, w = columns(trajectory(sys, 1000, Ttr = 10000));
test = SurrogateTest(RMCD(r = 0.5); rng)
α = 0.05
@test pvalue(independence(test, x, z, y)) < α
