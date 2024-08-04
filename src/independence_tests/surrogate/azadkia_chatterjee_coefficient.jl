function independence(test::SurrogateAssociationTest{<:AzadkiaChatterjeeCoefficient}, 
        x::AbstractVector, y_and_possibly_z...)
    (; est_or_measure, rng, surrogate, nshuffles) = test

    # Create a new instance of the measure with pre-allocated values. We can then 
    # re-use this struct to avoid excessive allocations. 
    measure = AzadkiaChatterjeeCoefficient(x, y_and_possibly_z...; 
        theiler = test.est_or_measure.theiler)

    Î = association(measure, x, y_and_possibly_z...)
    sx = surrogenerator(x, surrogate, rng)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        Îs[b] = association(est_or_measure, sx(), y_and_possibly_z...)
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateAssociationTestResult(1 + length(y_and_possibly_z), Î, Îs, p, nshuffles)
end
