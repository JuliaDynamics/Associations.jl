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