using Random: shuffle
using StatsBase: sample
using Setfield

# function LocalPermutationTest(measure::TransferEntropy, est::Nothing, args...; kwargs...)
#     txt = "A valid estimator must be provided as second argument to "*
#     "`LocalPermutationTest` when using the `TEShannon` measure.\n" *
#         "Do e.g. LocalPermutationTest(TEShannon(), FPVP())"
#     throw(ArgumentError(txt))
# end

function independence(test::LocalPermutationTest{<:MultivariateInformationMeasureEstimator{<:TransferEntropy}}, x::AbstractVector...)
    est_or_measure, nshuffles = deepcopy(test.est_or_measure), test.nshuffles
    if !(length(x) == 3) && est isa TransferEntropyEstimator
        msg = "`LocalPermutationTest` is not defined for pairwise transfer entropy with " *
        " `TransferEntropyEstimators`. " * 
            "Either provide a third timeseries to condition on, or use some other estimator."
        throw(ArgumentError(msg))
    end

    def = est_or_measure.definition
    # Below, the T variable also includes any conditional variables.
    S, T, T⁺, C = individual_marginals_te(def.embedding, x...)
    TC = StateSpaceSet(T, C)
    @assert length(T⁺) == length(S) == length(TC)
    N = length(x)

    X, Y = S, T⁺ # The source marginal `S` is the one being shuffled.
    Z = TC # The conditional variable
    est = convert_to_cmi_estimator(est_or_measure)
   
    Î = association(est, X, Y, Z)
    # works until here.

    Îs = permuted_Îs_te(S, T, T⁺, C, est, test)
    # TODO: make compatible with TransferEntropyEstimators.
    # THis requires a new permuted_Îs_te dedicated for that.
    # if est_or_measure isa TransferEntropyEstimator
    #     # @show "lll"
    #     # @show "heyo"
    #     # Î = association(est_or_measure, S, T, T⁺, C)
    #             #Îs = permuted_Îs_te(S, T, T⁺, C, est_or_measure, test)

    #     X, Y = S, T⁺ # The source marginal `S` is the one being shuffled.
    #     Z = TC # The conditional variable
    #     Îs = permuted_Îs(X, Y, Z, est_or_measure, test)

    # end

    p = count(Î .<= Îs) / nshuffles
    return LocalPermutationTestResult(length(x), Î, Îs, p, nshuffles)
end

# Runge's local permutation test can't be directly translated to transfer entropy specific 
# estimators like `Lindner`. However, but we can use a similar principle where 
# the source marginal `S` is shuffled according to local closeness in the 
# conditional marginal `C`. The `T` and `T⁺` marginals (i.e. all information)
# about the target variable is left untouched.
function permuted_Îs_te(S, T, T⁺, C, est_or_measure, test)
    rng, kperm, nshuffles, replace, w = test.rng, test.kperm, test.nshuffles, test.replace, test.w
    progress = ProgressMeter.Progress(nshuffles;
        desc = "LocalPermutationTest:",
        enabled = test.show_progress
    )
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
        Îs[n] = association(est_or_measure, Ŝ, T, T⁺, C)
        ProgressMeter.next!(progress)
    end
    return Îs
end
