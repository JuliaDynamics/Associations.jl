using StateSpaceSets
using Distances: Chebyshev

rng = Xoshiro(1234)
x = rand(rng, 200)
y = rand(rng, 200)
X = rand(rng, 200, 2) |> StateSpaceSet
Y = rand(rng, 200, 3) |> StateSpaceSet
@test_throws UndefKeywordError MCR()
@test association(MCR(; r = 0.5), x, y) isa Real
@test association(MCR(; r = 0.5), x, Y) isa Real
@test association(MCR(; r = 0.5), X, Y) isa Real
@test association(MCR(; r = 0.5, metric = Chebyshev()), x, y) isa Real
@test association(MCR(; r = 0.5, metric = Chebyshev()), X, Y) isa Real


test = SurrogateAssociationTest(MCR(r = 0.2); rng)
α = 0.05
# We should not be able to reject null for independent variables
@test pvalue(independence(test, x, y)) >= α
@test pvalue(independence(test, X, Y)) >= α
@test pvalue(independence(test, x, Y)) >= α

# For dependent variables, we should be able to reject the null hypothesis of
# independence. This goes both ways.
z = x .+ y
@test pvalue(independence(test, x, z)) < α
@test pvalue(independence(test, z, x)) < α

# Romano et al. claim that if A drives B, then ΔM(A | B) = M(A | B) - M(B | A) > 0
m = MCR(; r = 0.5)
Δxz = association(m, x, z) - association(m, z, x)
Δyz = association(m, y, z) - association(m, y, x)
@test Δxz > 0
@test Δyz > 0
