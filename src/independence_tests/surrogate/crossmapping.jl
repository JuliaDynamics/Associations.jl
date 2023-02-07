using Statistics: mean

function independence(test::SurrogateTest{<:CrossmapMeasure, <:CrossmapEstimator{Int}}, x, y)
    (; measure, est, rng, surrogate, nshuffles) = test
    Î = crossmap(measure, est, x, y)
    sx = surrogenerator(x, surrogate, rng)
    sy = surrogenerator(y, surrogate, rng)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        Îs[b] = crossmap(measure, est, sx(), sy())
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateTestResult(Î, Îs, p, nshuffles)
end

function independence(test::SurrogateTest{<:Ensemble{<:CrossmapMeasure, <:RandomVectors{Int}}}, x, y)
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
    return SurrogateTestResult(mean(Î), Îs, p, nshuffles)
end


# Independence tests are currently only defined for estimators operating on a single 
# library size.
const INVALID_ENSEMBLE = Ensemble{
    <:CrossmapMeasure, 
    <:CrossmapEstimator{<:Union{AbstractVector, AbstractRange}
    }}
const INVALID_CM_TEST = SurrogateTest{<:INVALID_ENSEMBLE}

function SurrogateTest(measure::CrossmapMeasure, est::CrossmapEstimator{<:Union{AbstractVector, AbstractRange}}, args...; kwargs...)
    T = typeof(est)
    txt = "\n`SurrogateTest` not implemented for estimator $T. Specifically,\n" *
        "`SurrogateTest(CCM(), RandomVectors(libsizes = 100:200:500, replace = true)))`" * 
        " will not work.\n" *
        "The estimator must operate on a single library size, e.g.\n" *
        "`SurrogateTest(CCM(), RandomVectors(libsizes = 100, replace = true))`.\n" 
       
    throw(ArgumentError(txt))
end

function SurrogateTest(e::INVALID_ENSEMBLE, args...; kwargs...)
    T = typeof(e.est)
    txt = "\n`SurrogateTest` not implemented for estimator $T. Specifically,\n" *
        "`SurrogateTest(CCM(), RandomVectors(libsizes = 100:200:500, replace = true)))`" * 
        " will not work.\n" *
        "The estimator must operate on a single library size, e.g.\n" *
        "`SurrogateTest(CCM(), RandomVectors(libsizes = 100, replace = true))`.\n" 
       
    throw(ArgumentError(txt))
end