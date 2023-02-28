using Random: shuffle
using StatsBase: sample

function LocalPermutationTest(measure::TransferEntropy, est::Nothing, args...; kwargs...)
    txt = "A valid estimator must be provided as second argument to "*
    "`LocalPermutationTest` when using the `TEShannon` measure.\n" *
        "Do e.g. LocalPermutationTest(TEShannon(), FPVP())"
    throw(ArgumentError(txt))
end

function independence(test::LocalPermutationTest{<:TransferEntropy{<:E}}, x::AbstractVector...) where E
    measure, est, nshuffles = test.measure, test.est, test.nshuffles
    # Below, the T variable also includes any conditional variables.
    S, T, T⁺, C = individual_marginals_te(measure.embedding, x...)
    TC = StateSpaceSet(T, C)
    @assert length(T⁺) == length(S) == length(TC)
    N = length(x)

    X, Y = T⁺, S
    Z = TC # The conditional variable
    cmi = te_to_cmi(measure)
    Î = estimate(cmi, est, X, Y, Z)
    Îs = permuted_Îs(X, Y, Z, cmi, est, test)
    p = count(Î .<= Îs) / nshuffles
    return LocalPermutationTestResult(length(x), Î, Îs, p, nshuffles)
end
