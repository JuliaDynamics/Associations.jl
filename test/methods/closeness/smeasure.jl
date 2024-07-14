x, y = rand(200), rand(200)
X, Y = StateSpaceSet(rand(200, 3)), StateSpaceSet(rand(200, 2))
Z, W = StateSpaceSet(rand(210, 2)), StateSpaceSet(rand(190, 4))
dx, τx = 2, 1
dy, τy = 2, 1

@test SMeasure() isa ClosenessMeasure
@test association(SMeasure(), x, y) isa Float64
@test association(SMeasure(dx = dx, τx = τx), x, Y) isa Float64
@test association(SMeasure(dy = dy, τy = τy), X, y) isa Float64
@test association(SMeasure(), X, Y) isa Float64
# test that multivariate StateSpaceSets are being length-matched
@test association(SMeasure(), X, Z) isa Float64
@test association(SMeasure(), W, X) isa Float64

# Deprecations
@test_logs (:warn, "Convenience function `s_measure` is deprecated. Use `association(SMeasure(; kwargs...), x, y)` instead.") s_measure(SMeasure(), x, y)
@test_logs (:warn, "Convenience function `s_measure` is deprecated. Use `association(SMeasure(; kwargs...), x, y)` instead.") s_measure(x, y)