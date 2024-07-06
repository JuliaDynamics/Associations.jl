# For `HLMS` measures (e.g. ``), mixtures of `Vector` and `StateSpaceSet`s are allowed,
# so we need to dispatch explicitly and call `s_measure` manually
# to avoid automatic conversion to `StateSpaceSet`s (which would ignore
# embedding parameters if input data are `Vector`s).
function independence(test::SurrogateAssociationTest{<:HLMS}, x, y)
    (; est_or_measure, rng, surrogate, nshuffles) = test
    Î = association(est_or_measure, x, y)
    sx = surrogenerator(x, surrogate, rng)
    sy = surrogenerator(y, surrogate, rng)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        Îs[b] = association(est_or_measure, sx(), sy())
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateAssociationTestResult(2, Î, Îs, p, nshuffles)
end
