using Random
rng = MersenneTwister(1234)
x, y = rand(rng, 500), rand(rng, 500)
z = x.+ y

d = 2
τ = -1

# Regular variant.
test_ccm = SurrogateTest(CCM(; d, τ), RandomVectors(libsizes = 300; replace = true))
test_pai = SurrogateTest(PAI(; d, τ), RandomVectors(libsizes = 300; replace = true))
@test_throws ArgumentError SurrogateTest(Ensemble(CCM(), RandomVectors(libsizes = 100:100:300)))
@test_throws ArgumentError SurrogateTest(CCM(), RandomVectors(libsizes = 100:100:300))

α = 0.03 # arbitrarily set confidence level to 1 - α
@test pvalue(independence(test_ccm, x, y)) > α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_pai, x, y)) > α

# Ensemble variant.
eccm = Ensemble(CCM(; d, τ), RandomVectors(libsizes = 100; replace = true))
epai = Ensemble(PAI(; d, τ), RandomVectors(libsizes = 100; replace = true))
test_ccm = SurrogateTest(eccm)
test_pai = SurrogateTest(epai)
@test pvalue(independence(test_ccm, x, y)) > α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_ccm, x, z)) < α
@test pvalue(independence(test_pai, x, y)) > α
