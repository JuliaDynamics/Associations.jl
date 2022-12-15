
using Distributions: MvNormal
using StateSpaceSets: Dataset
using Statistics: mean
using Random; rng = MersenneTwister(1234567)

################################################################
# Analytical tests.
#
# For same-dimensional multivariate Gaussians, the relative
# entropy is known analytically. We compute it using
# `divergence(Shannon(; base), Nx, Ny)` where `Nx` and `Ny`
# are `MvNormal`s.
################################################################

# With identical means covariance matrices, there should be zero relative entropy
σx = 0.5; σy = 0.5; vx = 1.0; vy = 1.0
Nx = MvNormal([0, 0], [vx σx; σx vx])
Ny = MvNormal([0, 0], [vy σy; σy vy])
Dx = Dataset([rand(rng, Nx) for i = 1:100000])
Dy = Dataset([rand(rng, Ny) for i = 1:100000])

re_true = CausalityTools.divergence(Shannon(), Nx, Ny)
re_wang = divergence(RelativeEntropyShannon(), Wang(k = 5, l = 5), Dx, Dy)
re_wang_trans = divergence(RelativeEntropyShannon(), WangTransformed(k = 5, l = 5), Dx, Dy)

# Just test that we're "resonably" close to zero.
@test -1e2 < round(re_wang, digits = 2) < 1e2
@test -1e2 < round(re_wang_trans, digits = 2) < 1e2

# Test another case with some randomly selected parameters
σx = 0.5
σy = 0.7
vx = 1.0
vy = 1.0
Nx = MvNormal([1, 1], [vx σx; σx vx])
Ny = MvNormal([0, 0], [vy σy; σy vy])
nreps = 50
re_wang = zeros(nreps)
re_wang_trans = zeros(nreps)
for i = 1:nreps
    local Dx = Dataset([rand(Nx) for i = 1:10000])
    local Dy = Dataset([rand(Ny) for i = 1:10000])
    re_wang[i] = divergence(RelativeEntropyShannon(), Wang(k = 5, l = 5), Dx, Dy)
    re_wang_trans[i] = divergence(RelativeEntropyShannon(), WangTransformed(k = 5, l = 5), Dx, Dy)
end

# Test that on average, we're within 10% of the target.
# Note: this does not at all hold for all distributions.
# The estimator is sensitive to the underlying distribution.
htrue = divergence(Shannon(), Nx, Ny)
@test htrue*0.9 < mean(re_wang) < htrue*1.1
@test htrue*0.9 < mean(re_wang_trans) < htrue*1.1
