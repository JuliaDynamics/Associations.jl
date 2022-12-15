using Distributions: MvNormal
using Random; rng = MersenneTwister(1234567)

# With identical means covariance matrices, there should be zero relative entropy
σx = 0.5; σy = 0.5; vx = 1.0; vy = 1.0
Nx = MvNormal([0, 0], [vx σx; σx vx])
Ny = MvNormal([0, 0], [vy σy; σy vy])
Dx = Dataset([rand(rng, Nx) for i = 1:100000])
Dy = Dataset([rand(rng, Ny) for i = 1:100000])

q = 1.5
base = 2
re_true = CausalityTools.divergence(Renyi(; q, base), Nx, Ny)
re_ps = divergence(RelativeEntropyRenyi(; q, base), PoczosSchneiderRE(k = 20), Dx, Dy)

# Just test that we're "resonably" close to zero
@test min(-0.1, re_true * 0.5) < re_ps < max(0.1, re_true * 1.5)
