x, y = rand(200), rand(200)
X, Y = StateSpaceSet(rand(200, 3)), StateSpaceSet(rand(200, 2))
Z, W = StateSpaceSet(rand(210, 2)), StateSpaceSet(rand(190, 4))
dx, τx = 2, 1
dy, τy = 2, 1

@test MMeasure() isa ClosenessMeasure
@test association(MMeasure(), x, y) isa Float64
@test association(MMeasure(dx = dx, τx = τx), x, Y) isa Float64
@test association(MMeasure(dy = dy, τy = τy), X, y) isa Float64
@test association(MMeasure(), X, Y) isa Float64
# test that multivariate StateSpaceSets are being length-matched
@test association(MMeasure(), X, Z) isa Float64
@test association(MMeasure(), W, X) isa Float64

# Deprecations
@test_logs (:warn, "Convenience function `m_measure` is deprecated. Use `association(MMeasure(; kwargs...), source, target) instead.") m_measure(MMeasure(), x, y)
@test_logs (:warn, "Convenience function `m_measure` is deprecated. Use `m_measure(MMeasure(; kwargs...), source, target)` instead.") m_measure(x, y)