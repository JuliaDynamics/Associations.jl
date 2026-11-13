using Test
using Associations
using Random
using LinearAlgebra
using Statistics

rng = Xoshiro(1234)

# ---------------------------------------------------------------------------------------
# Exactly balanced discrete systems, for which the plug-in estimates equal the true values.
# ---------------------------------------------------------------------------------------
# All 2^k combinations of k fair bits, each appearing exactly once.
function all_bit_combinations(k)
    return [[(i >> (j - 1)) & 1 for i in 0:(2^k - 1)] for j in 1:k]
end

# n-1 independent fair bits plus their parity.
function parity_system(n)
    bits = all_bit_combinations(n - 1)
    parity = reduce((a, b) -> a .⊻ b, bits)
    return (bits..., parity)
end

# n identical copies of one fair bit.
copy_system(n) = ntuple(_ -> [0, 1], n)

# n independent fair bits.
independent_system(n) = Tuple(all_bit_combinations(n))

disc = CodifyVariables(UniqueElements())
discrete_estimators(def) = [
    JointProbabilities(def, disc),
    EntropyDecomposition(def, PlugIn(Shannon()), disc),
]

# Values in bits: (TC, DTC, Ω, S).
# Parity of n bits:   TC = 1,     DTC = n - 1,  Ω = 2 - n,  S = n.
# n copies of a bit:  TC = n - 1, DTC = 1,      Ω = n - 2,  S = n.
# Independent bits:   all zero.
expected(::typeof(parity_system), n) = (1, n - 1, 2 - n, n)
expected(::typeof(copy_system), n) = (n - 1, 1, n - 2, n)
expected(::typeof(independent_system), n) = (0, 0, 0, 0)

@testset "Exact discrete systems" begin
    for system in (parity_system, copy_system, independent_system), n in 3:5
        x = system(n)
        tc, dtc, Ω, s = expected(system, n)
        for (def, val) in zip(
                (TotalCorrelation(), DualTotalCorrelation(), OInformation(), SInformation()),
                (tc, dtc, Ω, s))
            for est in discrete_estimators(def)
                @test association(est, x...) ≈ val atol = 1e-12
            end
        end
    end
end

@testset "Units" begin
    x = parity_system(3)
    @test association(JointProbabilities(OInformation(base = ℯ), disc), x...) ≈ -log(2)
    @test association(JointProbabilities(TotalCorrelation(base = ℯ), disc), x...) ≈ log(2)
end

@testset "Identities" begin
    x, y, z, w = rand(rng, 1:3, 2000), rand(rng, 1:4, 2000), rand(rng, 1:2, 2000), rand(rng, 1:3, 2000)
    z = z .+ (x .== y)
    vars = (x, y, z, w)
    tc = association(JointProbabilities(TotalCorrelation(), disc), vars...)
    dtc = association(JointProbabilities(DualTotalCorrelation(), disc), vars...)
    Ω = association(JointProbabilities(OInformation(), disc), vars...)
    s = association(JointProbabilities(SInformation(), disc), vars...)
    @test Ω ≈ tc - dtc
    @test s ≈ tc + dtc

    # Invariant under reordering of the variables.
    @test association(JointProbabilities(OInformation(), disc), w, z, x, y) ≈ Ω

    # For two variables, TC = DTC = mutual information and Ω = 0.
    mi = association(JointProbabilities(MIShannon(), disc), x, y)
    @test association(JointProbabilities(TotalCorrelation(), disc), x, y) ≈ mi
    @test association(JointProbabilities(DualTotalCorrelation(), disc), x, y) ≈ mi
    @test association(JointProbabilities(OInformation(), disc), x, y) ≈ 0 atol = 1e-12

    # For three variables, Ω equals the co-information I(X;Y) - I(X;Y|Z).
    cmi = association(JointProbabilities(CMIShannon(), disc), x, y, z)
    @test association(JointProbabilities(OInformation(), disc), x, y, z) ≈ mi - cmi

    # One input is not enough.
    @test_throws ArgumentError association(JointProbabilities(OInformation(), disc), x)
end

# ---------------------------------------------------------------------------------------
# Bias-corrected discrete estimators. Each subset entropy is estimated from the counts of
# that subset's observed joint outcomes; the references below are closed-form estimates
# from those counts, written independently of ComplexityMeasures.
# ---------------------------------------------------------------------------------------
using StatsBase: countmap

plugin_h(n) = (N = sum(n); -sum(k / N * log(k / N) for k in n))
miller_madow_h(n) = plugin_h(n) + (length(n) - 1) / (2 * sum(n))
horvitz_thompson_h(n) = (N = sum(n); -sum(k / N * log(k / N) / (1 - (1 - k / N)^N) for k in n))
# Chao & Shen (2003): Horvitz–Thompson with probabilities scaled by the coverage C = 1 - f₁/N.
function chao_shen_h(n)
    N = sum(n)
    C = 1 - count(==(1), n) / N
    return -sum(C * k / N * log(C * k / N) / (1 - (1 - C * k / N)^N) for k in n)
end
function jackknife_h(n)
    N = sum(n)
    loo = sum(n[j] * plugin_h(filter(>(0), [i == j ? k - 1 : k for (i, k) in enumerate(n)]))
        for j in eachindex(n))
    return N * plugin_h(n) - (N - 1) / N * loo
end
# Digamma at positive multiples of 1/2, by recurrence from ψ(1) and ψ(1/2).
function digamma_half(x)
    y, ψ = isinteger(x) ? (1.0, -Base.MathConstants.γ) : (0.5, -Base.MathConstants.γ - 2log(2))
    while y < x
        ψ += 1 / y
        y += 1
    end
    return ψ
end
# Grassberger (2003); equal to Schürmann (2004) with a = 1/2.
grassberger_G(k) = digamma_half(k) + 0.5 * (-1)^k * (digamma_half((k + 1) / 2) - digamma_half(k / 2))
grassberger_h(n) = (N = sum(n); digamma_half(N) - sum(k * grassberger_G(k) for k in n) / N)
dirichlet_h(n, a) = (q = (n .+ a) ./ (sum(n) + length(n) * a); -sum(q .* log.(q)))
# James–Stein shrinkage towards uniform over the observed outcomes (Hausser & Strimmer 2009).
function shrinkage_h(n)
    N, m = sum(n), length(n)
    p = n ./ N
    λ = clamp((1 - sum(abs2, p)) / ((N - 1) * sum(abs2, 1 / m .- p)), 0, 1)
    q = λ / m .+ (1 - λ) .* p
    return -sum(q .* log.(q))
end

@testset "Bias-corrected discrete estimators" begin
    rng_bc = Xoshiro(99)
    N = 120
    a, b = rand(rng_bc, 0:2, N), rand(rng_bc, 0:2, N)
    c = mod.(a .+ b .+ (rand(rng_bc, N) .< 0.2), 3)
    d = (a .== c) .+ rand(rng_bc, 0:1, N)
    vars = (a, b, c, d)
    M = hcat(vars...)
    subset_counts(S) = collect(values(countmap([Tuple(r) for r in eachrow(M[:, S])])))
    reference(h, def) = Associations.entropy_combination(def,
        Associations.SubsetEntropies(S -> h(subset_counts(S))), collect(1:4))

    cases = [
        (MillerMadow(Shannon()), RelativeAmount(), miller_madow_h),
        (HorvitzThompson(Shannon()), RelativeAmount(), horvitz_thompson_h),
        (ChaoShen(Shannon()), RelativeAmount(), chao_shen_h),
        (Jackknife(Shannon()), RelativeAmount(), jackknife_h),
        (GeneralizedSchuermann(Shannon(); a = 1.0), RelativeAmount(), grassberger_h),
        (Schuermann(Shannon(); a = 0.5), RelativeAmount(), grassberger_h),
        (PlugIn(Shannon()), BayesianRegularization(a = 0.5), n -> dirichlet_h(n, 0.5)),
    ]
    # The sample is small enough that some joint outcomes are singletons, so the Chao–Shen
    # coverage correction is exercised.
    @test any(==(1), subset_counts([1, 2, 3, 4]))
    for def in (TotalCorrelation(base = ℯ), DualTotalCorrelation(base = ℯ),
            OInformation(base = ℯ), SInformation(base = ℯ))
        for (est, pest, h) in cases
            @test association(EntropyDecomposition(def, est, disc, pest), vars...) ≈ reference(h, def)
        end
    end

    # The correction actually changes the estimate relative to plug-in.
    plug = association(EntropyDecomposition(OInformation(base = ℯ), PlugIn(Shannon()), disc), vars...)
    mm = association(EntropyDecomposition(OInformation(base = ℯ), MillerMadow(Shannon()), disc), vars...)
    @test !(plug ≈ mm)

    # Known bugs in ComplexityMeasures 3.12.1, recorded so they show up when fixed upstream.
    # 1. MillerMadow adds its (m - 1)/2N correction in nats whatever the base.
    @test_broken association(EntropyDecomposition(OInformation(), MillerMadow(Shannon()), disc), vars...) ≈
        association(EntropyDecomposition(OInformation(base = ℯ), MillerMadow(Shannon()), disc), vars...) / log(2)
    # 2. Shrinkage multiplies by Σ(tₖ - pₖ)² in λ instead of dividing by it.
    @test_broken association(EntropyDecomposition(OInformation(base = ℯ), PlugIn(Shannon()), disc, Shrinkage()), vars...) ≈
        reference(shrinkage_h, OInformation(base = ℯ))
    # 3. ChaoShen gives NaN when every observation is a singleton.
    @test_broken isfinite(association(EntropyDecomposition(OInformation(), ChaoShen(Shannon()), disc), 1:10, 1:10, 1:10))
end

# ---------------------------------------------------------------------------------------
# Gaussian systems with closed-form values.
# ---------------------------------------------------------------------------------------
# Gaussian entropy, up to constants that cancel in every measure here: ½ log det Σ_S.
gaussian_h(Σ, S) = 0.5 * logdet(Σ[S, S])
function gaussian_truth(Σ, groups)
    n = length(groups)
    H(gs) = gaussian_h(Σ, reduce(vcat, groups[gs]))
    hall = H(1:n)
    tc = sum(H([i]) for i in 1:n) - hall
    dtc = sum(H(setdiff(1:n, i)) for i in 1:n) - (n - 1) * hall
    return (tc, dtc, tc - dtc, tc + dtc) ./ log(2)
end

function sample_gaussian(rng, Σ, N)
    L = cholesky(Symmetric(Σ)).L
    return Matrix((L * randn(rng, size(Σ, 1), N))')
end

defs = (TotalCorrelation(), DualTotalCorrelation(), OInformation(), SInformation())

@testset "Gaussian estimator" begin
    A = randn(rng, 4, 4)
    Σ = A * A' + 0.5I
    N = 200_000
    X = sample_gaussian(rng, Σ, N)
    cols = collect(eachcol(X))
    truth = gaussian_truth(Σ, [[1], [2], [3], [4]])
    for (def, val) in zip(defs, truth)
        @test association(GaussianDecomposition(def), cols...) ≈ val atol = 0.02
        # The copula transform makes the estimate invariant to monotone transformations.
        warped = (exp.(cols[1]), cols[2] .^ 3, cols[3], -cols[4])
        @test association(GaussianDecomposition(def; copula = true), warped...) ≈ val atol = 0.03
    end

    # Grouped variables: the inputs are StateSpaceSets with several columns each.
    groups = [[1, 2], [3], [4]]
    truth_grouped = gaussian_truth(Σ, groups)
    G = (StateSpaceSet(X[:, 1:2]), X[:, 3], X[:, 4])
    for (def, val) in zip(defs, truth_grouped)
        @test association(GaussianDecomposition(def), G...) ≈ val atol = 0.02
    end

    # Two variables: TC equals the Gaussian mutual information.
    mi = association(GaussianMI(MIShannon(); normalize = false), StateSpaceSet(X[:, 1:2]), X[:, 3])
    @test association(GaussianDecomposition(TotalCorrelation()), StateSpaceSet(X[:, 1:2]), X[:, 3]) ≈ mi
end

@testset "k-NN entropy decomposition" begin
    A = randn(rng, 3, 3)
    Σ = A * A' + 0.5I
    # Each subset entropy is estimated separately, so the estimator's bias in 1, 2 and 3
    # dimensions does not cancel; at N = 20 000 the error is a few hundredths of a bit.
    X = sample_gaussian(rng, Σ, 20_000)
    truth = gaussian_truth(Σ, [[1], [2], [3]])
    for (def, val) in zip(defs, truth)
        est = EntropyDecomposition(def, Kraskov(Shannon(); k = 5))
        @test association(est, eachcol(X)...) ≈ val atol = 0.08
    end
end

@testset "Binned continuous data" begin
    A = randn(rng, 3, 3)
    Σ = A * A' + 0.5I
    X = sample_gaussian(rng, Σ, 200_000)
    c = collect(eachcol(X))
    truth = gaussian_truth(Σ, [[1], [2], [3]])
    d = CodifyVariables(ValueBinning(RectangularBinning(16)))
    for (def, val) in zip(defs, truth)
        joint = association(JointProbabilities(def, d), c...)
        separate = association(EntropyDecomposition(def, PlugIn(Shannon()), d), c...)
        # The two discrete routes bin identically, so they agree exactly.
        @test joint ≈ separate
        # With enough samples per bin, binning recovers the continuous value.
        @test joint ≈ val atol = 0.06
    end
end

@testset "Grouped inputs with ValueBinning" begin
    # Regression test: `codified_marginals` with `ValueBinning` used to give each input a
    # single column of the joint binning, so a grouped input silently lost its other
    # columns and the inputs after it were given the wrong columns.
    N = 20_000
    a, b, c = randn(rng, N), randn(rng, N), randn(rng, N)
    y = a .+ 0.1 .* randn(rng, N) # depends on the first column of the group only
    G = StateSpaceSet(a, b)
    d = CodifyVariables(ValueBinning(RectangularBinning(4)))

    enc = codified_marginals(d, G, y)
    @test length(enc) == 2
    @test dimension(StateSpaceSet(enc[1])) == 2

    for def in defs
        joint = association(JointProbabilities(def, d), G, y, c)
        separate = association(EntropyDecomposition(def, PlugIn(Shannon()), d), G, y, c)
        @test joint ≈ separate
    end
    # Column `a` of the group carries the dependence, so it must be seen.
    @test association(EntropyDecomposition(TotalCorrelation(), PlugIn(Shannon()), d), G, y) > 0.5

    # The same encodings feed the discrete mutual information.
    @test association(EntropyDecomposition(MIShannon(), PlugIn(Shannon()), d), G, y) ≈
        association(JointProbabilities(MIShannon(), d), G, y)
end

# ---------------------------------------------------------------------------------------
# Gradients of the O-information.
# ---------------------------------------------------------------------------------------
@testset "O-information gradients" begin
    # XOR triplet (x, y, z) plus an independent bit w, all 16 combinations once.
    bits = all_bit_combinations(3)
    x, y, w = bits
    z = x .⊻ y
    est = JointProbabilities(OInformation(), disc)

    g1 = oinformation_gradients(est, x, y, z, w)
    @test length(g1) == 4
    @test g1[1:3] ≈ fill(-1.0, 3) atol = 1e-12  # removing any XOR member destroys the synergy
    @test g1[4] ≈ 0 atol = 1e-12                # the independent bit contributes nothing

    g2 = oinformation_gradients(est, x, y, z, w; order = 2)
    @test size(g2) == (4, 4)
    @test g2 ≈ g2'
    # ∂²_ij Ω = Ω(X) - Ω(X₋ᵢ) - Ω(X₋ⱼ) + Ω(X₋ᵢⱼ)
    Ω(vs...) = length(vs) < 3 ? 0.0 : association(est, vs...)
    @test g2[1, 2] ≈ Ω(x, y, z, w) - Ω(y, z, w) - Ω(x, z, w) + Ω(z, w) atol = 1e-12

    # The Gaussian estimator's gradients agree with direct subset evaluation.
    A = randn(rng, 4, 4)
    X = sample_gaussian(rng, A * A' + 0.5I, 2000)
    c = collect(eachcol(X))
    gest = GaussianDecomposition(OInformation())
    gg = oinformation_gradients(gest, c...)
    @test gg[2] ≈ association(gest, c...) - association(gest, c[1], c[3], c[4])

    @test_throws ArgumentError oinformation_gradients(est, x, y, z; order = 3)
    @test_throws ArgumentError oinformation_gradients(JointProbabilities(TotalCorrelation(), disc), x, y, z)
end
