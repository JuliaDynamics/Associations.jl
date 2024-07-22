using Test
using Random
rng = MersenneTwister(1234)
n = 250
x, y = rand(rng, n), rand(rng, n)
z = x.+ y

d = 2
τ = -1

# Regular variant.
est_ccm = RandomVectors(CCM(; d, τ); libsizes = 100, replace = true, rng)
test_ccm = SurrogateAssociationTest(est_ccm; rng, nshuffles = 19)

est_pai = RandomVectors(PAI(; d, τ); libsizes = 100, replace = true, rng)
test_pai = SurrogateAssociationTest(est_pai; rng, nshuffles = 19)

# Invalid syntax.
@test_throws ArgumentError SurrogateAssociationTest(CCM(), RandomVectors(libsizes = 100:100:300))

α = 0.03 # arbitrarily set confidence level to 1 - α
@test pvalue(independence(test_ccm, x, y)) > α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_pai, x, y)) > α

# Ensemble variant.
eccm = Ensemble(RandomVectors(CCM(; d, τ); libsizes = 50, replace = true, rng))
epai = Ensemble(RandomVectors(PAI(; d, τ); libsizes = 50, replace = true, rng))
test_ccm = SurrogateAssociationTest(eccm; rng, nshuffles = 19)
test_pai = SurrogateAssociationTest(epai; rng, nshuffles = 19)
@test pvalue(independence(test_ccm, x, y)) > α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_pai, x, y)) > α
