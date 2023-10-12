using Statistics: mean

function independence(test::SurrogateAssociationTest{<:CrossmapMeasure, <:CrossmapEstimator{Int}}, x, y)
    (; measure, est, rng, surrogate, nshuffles) = test
    Î = crossmap(measure, est, x, y)
    sx = surrogenerator(x, surrogate, rng)
    sy = surrogenerator(y, surrogate, rng)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        Îs[b] = crossmap(measure, est, sx(), sy())
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateAssociationTestResult(2, Î, Îs, p, nshuffles)
end

function independence(test::SurrogateAssociationTest{<:Ensemble{<:CrossmapMeasure, <:RandomVectors{Int}}}, x, y)
    (; measure, est, rng, surrogate, nshuffles) = test
    Î = crossmap(measure, x, y) # A vector of length `measure.nreps`
    sx = surrogenerator(x, surrogate, rng)
    sy = surrogenerator(y, surrogate, rng)
    Îs = Vector{eltype(1.0)}(undef, 0)
    sizehint!(Îs, nshuffles * measure.nreps)

    for b in 1:nshuffles
        append!(Îs, crossmap(measure, sx(), sy()))
    end
    p = count(mean(Î) .<= Îs) / (nshuffles * measure.nreps)
    return SurrogateAssociationTestResult(2, mean(Î), Îs, p, nshuffles)
end


# Independence tests are currently only defined for estimators operating on a single
# library size.
const INVALID_ENSEMBLE = Ensemble{
    <:CrossmapMeasure,
    <:CrossmapEstimator{<:Union{AbstractVector, AbstractRange}
    }}
const INVALID_CM_TEST = SurrogateAssociationTest{<:INVALID_ENSEMBLE}

function SurrogateAssociationTest(measure::CrossmapMeasure, est::CrossmapEstimator{<:Union{AbstractVector, AbstractRange}}, args...; kwargs...)
    T = typeof(est)
    txt = "\n`SurrogateAssociationTest` not implemented for estimator $T. Specifically,\n" *
        "`SurrogateAssociationTest(CCM(), RandomVectors(libsizes = 100:200:500, replace = true)))`" *
        " will not work.\n" *
        "The estimator must operate on a single library size, e.g.\n" *
        "`SurrogateAssociationTest(CCM(), RandomVectors(libsizes = 100, replace = true))`.\n"

    throw(ArgumentError(txt))
end

function SurrogateAssociationTest(e::INVALID_ENSEMBLE, args...; kwargs...)
    T = typeof(e.est)
    txt = "\n`SurrogateAssociationTest` not implemented for estimator $T. Specifically,\n" *
        "`SurrogateAssociationTest(CCM(), RandomVectors(libsizes = 100:200:500, replace = true)))`" *
        " will not work.\n" *
        "The estimator must operate on a single library size, e.g.\n" *
        "`SurrogateAssociationTest(CCM(), RandomVectors(libsizes = 100, replace = true))`.\n"

    throw(ArgumentError(txt))
end
