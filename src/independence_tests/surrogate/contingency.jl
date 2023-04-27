
# Explicit dispatch for independence for `Contingency` estimator, because we don't
# want to convert categorical input data to `StateSpaceSets`.
function independence(test::SurrogateTest{MEASURE, <:Contingency}, x, y) where MEASURE
    (; measure, est, rng, surrogate, nshuffles) = test
    @assert length(x) == length(y)
    N = length(x)
    Î = estimate(measure, est, x, y)
    sx = surrogenerator(x, surrogate, rng)
    sy = surrogenerator(y, surrogate, rng)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        Îs[b] = estimate(measure, est, sx(), sy())
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateTestResult(2, Î, Îs, p, nshuffles)
end

function independence(test::SurrogateTest{MEASURE, <:Contingency}, x, y, z) where MEASURE
    (; measure, est, rng, surrogate, nshuffles) = test
    @assert length(x) == length(y) == length(z)
    N = length(x)
    Î = estimate(measure, est, x, y, z)
    s = surrogenerator(x, surrogate, rng)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        Îs[b] = estimate(measure, est, s(), y, z)
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateTestResult(3, Î, Îs, p, nshuffles)
end
