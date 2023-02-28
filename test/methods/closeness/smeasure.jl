x, y = rand(100), rand(100)
X, Y = StateSpaceSet(rand(100, 3)), StateSpaceSet(rand(100, 2))
Z, W = StateSpaceSet(rand(110, 2)), StateSpaceSet(rand(90, 4))
dx, τx = 2, 1
dy, τy = 2, 1

# Compat
@test s_measure(x, y) isa Float64
@test s_measure(x, Y, dx = dx, τx = τx) isa Float64
@test s_measure(X, y, dy = dy, τy = τy) isa Float64
@test s_measure(X, Y) isa Float64
# test that multivariate StateSpaceSets are being length-matched
@test s_measure(X, Z) isa Float64
@test s_measure(W, X) isa Float64

# V2.X
@test s_measure(SMeasure(), x, y) isa Float64
@test s_measure(SMeasure(dx = dx, τx = τx), x, Y) isa Float64
@test s_measure(SMeasure(dy = dy, τy = τy), X, y) isa Float64
@test s_measure(SMeasure(), X, Y) isa Float64
# test that multivariate StateSpaceSets are being length-matched
@test s_measure(SMeasure(), X, Z) isa Float64
@test s_measure(SMeasure(), W, X) isa Float64
