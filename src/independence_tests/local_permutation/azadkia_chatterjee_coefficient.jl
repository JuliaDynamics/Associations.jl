function independence(test::LocalPermutationTest{<:AzadkiaChatterjeeCoefficient}, x::AbstractVector, y, z)
    est_or_measure, nshuffles = test.est_or_measure, test.nshuffles
    # Make sure that the measure is compatible with the input data.
    verify_number_of_inputs_vars(est_or_measure, 3)

    Y, Z = StateSpaceSet(y), StateSpaceSet(z)
    @assert length(x) == length(Y) == length(Z)
    Î = association(est_or_measure, x, Y, Z)
    Îs = permuted_Îs(x, Y, Z, est_or_measure, test)
    p = count(Î .<= Îs) / nshuffles
    return LocalPermutationTestResult(3, Î, Îs, p, nshuffles)
end