using Random: shuffle
using StatsBase: sample

function LocalPermutationTest(measure::TransferEntropy, est::Nothing, args...; kwargs...)
    txt = "A valid estimator must be provided as second argument to "*
    "`LocalPermutationTest` when using the `TEShannon` measure.\n" *
        "Do e.g. LocalPermutationTest(TEShannon(), FPVP())"
    throw(ArgumentError(txt))
end

function independence(test::LocalPermutationTest{<:TransferEntropy{<:E}}, x::AbstractVector...) where E
    if !(length(x) == 3)
        msg = "`LocalPermutationTest` is not defined for pairwise transfer entropy. " * 
            "Two input timeseries were provided. Please provide a third timeseries to condition on."
        throw(ArgumentError(msg))
    end
    measure, est, nshuffles = test.measure, test.est, test.nshuffles
    # Below, the T variable also includes any conditional variables.
    S, T, T⁺, C = individual_marginals_te(measure.embedding, x...)
    TC = StateSpaceSet(T, C)
    @assert length(T⁺) == length(S) == length(TC)
    N = length(x)

    if est isa TransferEntropyEstimator
        Î = estimate(measure, est, S, T, T⁺, C)
        Îs = permuted_Îs_te(S, T, T⁺, C, measure, est, test)
    else
        X, Y = S, T⁺ # The source marginal `S` is the one being shuffled.
        Z = TC # The conditional variable
        cmi = te_to_cmi(measure)
        Î = estimate(cmi, est, X, Y, Z)
        Îs = permuted_Îs(X, Y, Z, cmi, est, test)
    end

    p = count(Î .<= Îs) / nshuffles
    return LocalPermutationTestResult(length(x), Î, Îs, p, nshuffles)
end

# Runge's local permutation test can't be directly translated to transfer entropy specific 
# estimators like `Lindner`. However, but we can use a similar principle where 
# the source marginal `S` is shuffled according to local closeness in the 
# conditional marginal `C`. The `T` and `T⁺` marginals (i.e. all information)
# about the target variable is left untouched.
function permuted_Îs_te(S, T, T⁺, C, measure::TransferEntropy, est, test)
    rng, kperm, nshuffles, replace, w = test.rng, test.kperm, test.nshuffles, test.replace, test.w

    N = length(S)
    test.kperm < N || throw(ArgumentError("kperm must be smaller than input data length"))

    # Search for neighbors in the conditional marginal
    tree_C = KDTree(C, Chebyshev())
    idxs_C = bulkisearch(tree_C, C, NeighborNumber(kperm), Theiler(w))

    # Shuffle source marginal `S` based on local closeness in C.
    Ŝ = deepcopy(S)
    Nᵢ = MVector{kperm, Int}(zeros(kperm)) # A statically sized copy
    πs = shuffle(rng, 1:N)
    Îs = zeros(nshuffles)
    for n in 1:nshuffles
        if replace
            shuffle_with_replacement!(Ŝ, S, idxs_C, rng)
        else
            shuffle_without_replacement!(Ŝ, S, idxs_C, kperm, rng, Nᵢ, πs)
        end
        Îs[n] = estimate(measure, est, Ŝ, T, T⁺, C)
    end
    return Îs
end
