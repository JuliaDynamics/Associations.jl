export Contingency

"""
    Contingency <: ProbabilitiesEstimator
    Contingency(est::Union{ProbabilitiesEstimator, Nothing} = nothing)

`Contingency` is a probabilities estimator that transforms input data to a multidimensional
probability mass function (internally represented as [`ContingencyMatrix`](@ref).

It works directly on raw discrete/categorical data. Alternatively, if a
[`ProbabilitiesEstimator`](@ref) `est` for which [`marginal_encodings`](@ref) is implemented
is given, then input data are first discretized before creating the contingency matrix.
"""
Base.@kwdef struct Contingency{E <: Union{Nothing, ProbabilitiesEstimator}} <: ProbabilitiesEstimator
    est::E = nothing
end

# Explicit dispatch for independence for `Contingency` estimator, because we don't
# want to convert categorical input data to `Datasets`.
function independence(test::SurrogateTest{MEASURE, <:Contingency}, x, y) where MEASURE
    (; measure, est, rng, surrogate, nsurr) = test
    @assert length(x) == length(y)
    N = length(x)
    Î = estimate(measure, est, x, y)
    sx = surrogenerator(x, surrogate, rng)
    sy = surrogenerator(y, surrogate, rng)
    Îs = zeros(nsurr)
    for b in 1:nsurr
        Îs[b] = estimate(measure, est, sx(), sy())
    end
    p = count(Î .<= Îs) / nsurr

    return SurrogateTestResult(Î, Îs, p, nsurr)
end

function independence(test::SurrogateTest{MEASURE, <:Contingency}, x, y, z) where MEASURE
    @show "heyo"
    (; measure, est, rng, surrogate, nsurr) = test
    @assert length(x) == length(y) == length(z)
    N = length(x)
    Î = estimate(measure, est, x, y, z)
    s = surrogenerator(x, surrogate, rng)
    Îs = zeros(nsurr)
    for b in 1:nsurr
        Îs[b] = estimate(measure, est, s(), y, z)
    end
    p = count(Î .<= Îs) / nsurr

    return SurrogateTestResult(Î, Îs, p, nsurr)
end
