x, y = rand(200), rand(200)
X, Y = StateSpaceSet(rand(200, 3)), StateSpaceSet(rand(200, 2))
Z, W = StateSpaceSet(rand(210, 2)), StateSpaceSet(rand(190, 4))
dx, τx = 2, 1
dy, τy = 2, 1

@test HMeasure() isa ClosenessMeasure
@test association(HMeasure(), x, y) isa Float64
@test association(HMeasure(dx = dx, τx = τx), x, Y) isa Float64
@test association(HMeasure(dy = dy, τy = τy), X, y) isa Float64
@test association(HMeasure(), X, Y) isa Float64
# test that multivariate StateSpaceSets are being length-matched
@test association(HMeasure(), X, Z) isa Float64
@test association(HMeasure(), W, X) isa Float64

# Deprecations
@test_logs (:warn, "Convenience function `h_measure` is deprecated. Use `association(HMeasure(; kwargs...), source, target) instead.") h_measure(HMeasure(), x, y)
@test_logs (:warn, "Convenience function `h_measure` is deprecated. Use `h_measure(HMeasure(; kwargs...), source, target)` instead.") h_measure(x, y)