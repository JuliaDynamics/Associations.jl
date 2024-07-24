using Random: shuffle
using StatsBase: sample

# ----------------------------------------------------------------
# For dedicated estimators of transfer entropy
# ----------------------------------------------------------------
function independence(test::LocalPermutationTest{<:TransferEntropyEstimator{<:TransferEntropy}}, x::AbstractVector...)
    est, nshuffles = deepcopy(test.est_or_measure), test.nshuffles
    if !(length(x) == 3) && est isa TransferEntropyEstimator
        msg = "`LocalPermutationTest` with $(typeof(est).name.name) is undefined for $(length(x)) " *
        "input time series. Please provide a third timeseries to condition on."
        throw(ArgumentError(msg))
    end

    def = est.definition
    # The source marginal `S` is the one being shuffled.
    S, T, T⁺, C = individual_marginals_te(def.embedding, x...)
    @assert length(T⁺) == length(S) == length(C) == length(T)
    Î = estimate_from_marginals(est, S, T, T⁺, C)
    Îs = permuted_Îs(S, T, T⁺, C, est, test)
    p = count(Î .<= Îs) / nshuffles
    return LocalPermutationTestResult(length(x), Î, Îs, p, nshuffles)
end
function permuted_Îs(S, T, T⁺, C, est::TransferEntropyEstimator, test::LocalPermutationTest)
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
        Îs[n] = estimate_from_marginals(est, Ŝ, T, T⁺, C)
        ProgressMeter.next!(progress)
    end
    return Îs
end

# ----------------------------------------------------------------
# For other estimators of transfer entropy
# ----------------------------------------------------------------
function independence(test::LocalPermutationTest{<:MultivariateInformationMeasureEstimator{<:TransferEntropy}}, x::AbstractVector...)
    est, nshuffles = deepcopy(test.est_or_measure), test.nshuffles
    if !(length(x) == 3)
        msg = "`LocalPermutationTest` with estimator $(typeof(est).name.name) is undefined for $(length(x)) " *
        "input time series. Please provide a third timeseries to condition on."
        throw(ArgumentError(msg))
    end
    # Below, the T variable also includes any conditional variables.
    S, T, T⁺, C = individual_marginals_te(est.definition.embedding, x...)
    TC = StateSpaceSet(T, C)
    @assert length(T⁺) == length(S) == length(TC)
    N = length(x)

    X, Y = S, T⁺ # The source marginal `S` is the one being shuffled.
    Z = TC # The conditional variable
    cmi_est = convert_to_cmi_estimator(est)
   
    Î = association(cmi_est, X, Y, Z)
    Îs = permuted_Îs_other(S, T, T⁺, C, cmi_est, test)
    p = count(Î .<= Îs) / nshuffles
    return LocalPermutationTestResult(length(x), Î, Îs, p, nshuffles)
end

# Runge's local permutation test can't be directly translated to transfer entropy specific 
# estimators like `Lindner`. However, but we can use a similar principle where 
# the source marginal `S` is shuffled according to local closeness in the 
# conditional marginal `C`. The `T` and `T⁺` marginals (i.e. all information)
# about the target variable is left untouched.
function permuted_Îs_other(S, T, T⁺, C, cmi_est, test)
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

    TC = StateSpaceSet(T, C)

    for n in 1:nshuffles
        if replace
            shuffle_with_replacement!(Ŝ, S, idxs_C, rng)
        else
            shuffle_without_replacement!(Ŝ, S, idxs_C, kperm, rng, Nᵢ, πs)
        end
        
        Îs[n] = association(cmi_est, Ŝ, T⁺, TC)
        ProgressMeter.next!(progress)
    end
    return Îs
end
