x, y = rand(200), rand(200)
X, Y = StateSpaceSet(rand(200, 3)), StateSpaceSet(rand(200, 2))
Z, W = StateSpaceSet(rand(210, 2)), StateSpaceSet(rand(190, 4))
dx, τx = 2, 1
dy, τy = 2, 1

@test LMeasure() isa ClosenessMeasure
@test association(LMeasure(), x, y) isa Float64
@test association(LMeasure(dx = dx, τx = τx), x, Y) isa Float64
@test association(LMeasure(dy = dy, τy = τy), X, y) isa Float64
@test association(LMeasure(), X, Y) isa Float64
# test that multivariate StateSpaceSets are being length-matched
@test association(LMeasure(), X, Z) isa Float64
@test association(LMeasure(), W, X) isa Float64

# Deprecations
@test_logs (:warn, "Convenience function `l_measure` is deprecated. Use `association(LMeasure(; kwargs...), source, target) instead.") l_measure(LMeasure(), x, y)
@test_logs (:warn, "Convenience function `l_measure` is deprecated. Use `l_measure(LMeasure(; kwargs...), source, target)` instead.") l_measure(x, y)