# For `SMeasure`, mixtures of `Vector` and `Dataset`s are allowed,
# so we need to dispatch explicitly and call `s_measure` manually
# to avoid automatic conversion to `Dataset`s (which would ignore
# embedding parameters if input data are `Vector`s).
function independence(test::SurrogateTest{<:SMeasure}, x, y)
    (; measure, est, rng, surrogate, nshuffles) = test
    Î = s_measure(measure, x, y)
    sx = surrogenerator(x, surrogate, rng)
    sy = surrogenerator(y, surrogate, rng)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        Îs[b] = s_measure(measure, sx(), sy())
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateTestResult(Î, Îs, p, nshuffles)
end