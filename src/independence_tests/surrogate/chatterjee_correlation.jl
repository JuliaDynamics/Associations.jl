function independence(test::SurrogateAssociationTest{<:ChatterjeeCorrelation}, 
        x::AbstractVector, y::AbstractVector)
    (; est_or_measure, rng, surrogate, nshuffles) = test

    # Create a new instance of the measure with pre-allocated values. We can then 
    # re-use this struct to avoid excessive allocations. 
    m = est_or_measure # the old measure
    measure = ChatterjeeCorrelation(x, y, rng = m.rng, lt = m.lt, handle_ties = m.handle_ties)

    Î = association(measure, x, y)
    sx = surrogenerator(x, surrogate, rng)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        Îs[b] = association(est_or_measure, sx(), y)
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateAssociationTestResult(2, Î, Îs, p, nshuffles)
end
