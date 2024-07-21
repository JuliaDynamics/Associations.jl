using Test
using Random
rng = MersenneTwister(1234)
x, y = rand(rng, 500), rand(rng, 500)
z = x.+ y

d = 2
τ = -1

# Regular variant.
test_ccm = SurrogateAssociationTest(RandomVectors(CCM(; d, τ); libsizes = 300, replace = true, rng))
test_pai = SurrogateAssociationTest(RandomVectors(PAI(; d, τ); libsizes = 300, replace = true, rng))
@test_throws ArgumentError SurrogateAssociationTest(Ensemble(CCM(), RandomVectors(libsizes = 100:100:300)))
@test_throws ArgumentError SurrogateAssociationTest(CCM(), RandomVectors(libsizes = 100:100:300))

α = 0.03 # arbitrarily set confidence level to 1 - α
@test pvalue(independence(test_ccm, x, y)) > α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_pai, x, y)) > α

# Ensemble variant.
eccm = Ensemble(RandomVectors(CCM(; d, τ); libsizes = 100, replace = true, rng))
epai = Ensemble(RandomVectors(PAI(; d, τ); libsizes = 100, replace = true, rng))
test_ccm = SurrogateAssociationTest(eccm)
test_pai = SurrogateAssociationTest(epai)
@test pvalue(independence(test_ccm, x, y)) > α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_pai, x, y)) > α
