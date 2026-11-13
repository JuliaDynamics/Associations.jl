using Test
using Associations
using Random
using LinearAlgebra
using Statistics

rng = Xoshiro(2718)
disc = CodifyVariables(UniqueElements())
pid_estimator(r; base = 2) = JointProbabilities(PartialInformationDecomposition(r; base), disc)
two_source_atoms(res) = (redundancy(res), unique_information(res, 1), unique_information(res, 2), synergy(res))

# ---------------------------------------------------------------------------------------
# The redundancy lattice.
# ---------------------------------------------------------------------------------------
@testset "Redundancy lattice" begin
    # The number of nodes is the Dedekind number minus 2.
    for (n, m) in zip(1:5, (1, 4, 18, 166, 7579))
        @test length(Associations.redundancy_lattice(n).nodes) == m
    end
    L = Associations.redundancy_lattice(3)
    m = length(L.nodes)
    bottom = L.index[sort(UInt[1, 2, 4])]  # {1}{2}{3}
    top = L.index[UInt[7]]                 # {123}
    @test all(L.below[bottom, j] for j in 1:m if j != bottom)
    @test all(L.below[i, top] for i in 1:m if i != top)
    @test first(L.order) == bottom && last(L.order) == top
    # The Möbius inversion reproduces any redundancy function from the atoms.
    r = randn(rng, m)
    a = Associations.mobius_inversion(L, r)
    @test all(r[j] ≈ a[j] + sum(a[i] for i in 1:m if L.below[i, j]; init = 0.0) for j in 1:m)
end

# ---------------------------------------------------------------------------------------
# Exact discrete cases. Every input combination appears equally often, so the plug-in
# estimate is the exact value. Values in bits: (redundancy, unique 1, unique 2, synergy).
# ---------------------------------------------------------------------------------------
@testset "Classic two-source gates" begin
    x1, x2 = [0, 0, 1, 1], [0, 1, 0, 1]
    hand = 0.5 * log2(4 / 3) + 0.25 * log2(1 / 3) + 0.25 * 2  # I(T; X1) for AND = 0.3113
    cases = [
        ("XOR", x1, x2, x1 .⊻ x2, (0, 0, 0, 1), (0, 0, 0, 1)),
        ("AND", x1, x2, x1 .& x2, (hand, 0, 0, 0.5), (hand, 0, 0, 0.5)),
        ("RDN", [0, 1], [0, 1], [0, 1], (1, 0, 0, 0), (1, 0, 0, 0)),
        ("UNQ", x1, x2, x1, (0, 1, 0, 0), (0, 1, 0, 0)),
        # Two independent bits copied into the target: both measures give 1 bit of
        # redundancy and 1 of synergy, the known shortcoming of I_min and I_MMI.
        ("COPY", x1, x2, 2 .* x1 .+ x2, (1, 0, 0, 1), (1, 0, 0, 1)),
    ]
    for (name, s1, s2, t, wb, mmi) in cases
        @test all(isapprox.(two_source_atoms(association(pid_estimator(IMin()), s1, s2, t)), wb; atol = 1e-12))
        @test all(isapprox.(two_source_atoms(association(pid_estimator(IMMI()), s1, s2, t)), mmi; atol = 1e-12))
    end
end

@testset "I_min and I_MMI differ (reference: dit)" begin
    counts = Dict((0, 0, 0) => 4, (0, 0, 1) => 1, (0, 1, 0) => 1, (0, 1, 2) => 3, (1, 0, 1) => 4,
        (1, 0, 0) => 1, (1, 1, 2) => 2, (1, 1, 1) => 2, (0, 1, 1) => 1)
    rows = reduce(vcat, [fill(o, c) for (o, c) in counts])
    s1, s2, t = getindex.(rows, 1), getindex.(rows, 2), getindex.(rows, 3)
    # dit 2.3: PID_WB and PID_MMI.
    wb = (0.111349894088, 0.084276647793, 0.279513802769, 0.131026674062)
    mmi = (0.195626541881, 0.0, 0.195237154976, 0.215303321854)
    @test all(isapprox.(two_source_atoms(association(pid_estimator(IMin()), s1, s2, t)), wb; atol = 1e-10))
    @test all(isapprox.(two_source_atoms(association(pid_estimator(IMMI()), s1, s2, t)), mmi; atol = 1e-10))
end

@testset "Three sources" begin
    bits = [[(i >> (j - 1)) & 1 for i in 0:7] for j in 1:3]
    parity = bits[1] .⊻ bits[2] .⊻ bits[3]
    for r in (IMin(), IMMI())
        res = association(pid_estimator(r), bits..., parity)
        @test synergy(res) ≈ 1 atol = 1e-12
        @test all(abs(v) < 1e-12 for (node, v) in atoms(res) if node != [[1, 2, 3]])
        # Three copies of one bit: everything is redundant.
        res = association(pid_estimator(r), [0, 1], [0, 1], [0, 1], [0, 1])
        @test redundancy(res) ≈ 1 atol = 1e-12
        @test sum(last, atoms(res)) ≈ 1 atol = 1e-12
        @test length(atoms(res)) == 18
    end
end

@testset "Identities" begin
    N = 3000
    x1, x2, x3 = rand(rng, 0:2, N), rand(rng, 0:1, N), rand(rng, 0:2, N)
    t = mod.(x1 .+ x2 .* x3 .+ (rand(rng, N) .< 0.2), 3)
    mi(xs...) = association(JointProbabilities(MIShannon(), disc), StateSpaceSet(xs...), t)
    for r in (IMin(), IMMI())
        res = association(pid_estimator(r), x1, x2, x3, t)
        # The atoms sum to the joint mutual information.
        @test sum(last, atoms(res)) ≈ mi(x1, x2, x3)
        # The atoms below a single source sum to that source's mutual information.
        res2 = association(pid_estimator(r), x1, x2, t)
        R, U1, U2, S = two_source_atoms(res2)
        @test R + U1 ≈ mi(x1)
        @test R + U2 ≈ mi(x2)
        @test R + U1 + U2 + S ≈ mi(x1, x2)
        # Swapping the sources swaps the unique informations.
        swapped = two_source_atoms(association(pid_estimator(r), x2, x1, t))
        @test all(swapped .≈ (R, U2, U1, S))
    end
    # I_MMI redundancy is the smallest single-source mutual information.
    res = association(pid_estimator(IMMI()), x1, x2, t)
    @test redundancy(res) ≈ min(mi(x1), mi(x2))
    # Units.
    x, y = [0, 0, 1, 1], [0, 1, 0, 1]
    @test synergy(association(pid_estimator(IMin(); base = ℯ), x, y, x .⊻ y)) ≈ log(2)
end

# ---------------------------------------------------------------------------------------
# Continuous data with I_MMI.
# ---------------------------------------------------------------------------------------
@testset "Gaussian I_MMI" begin
    N = 100_000
    Σx = [1.0 0.5; 0.5 1.0]
    X = Matrix((cholesky(Σx).L * randn(rng, 2, N))')
    β, σ2 = [0.8, 0.3], 0.36
    t = X * β .+ sqrt(σ2) .* randn(rng, N)
    x1, x2 = X[:, 1], X[:, 2]
    # Analytic decomposition, in bits.
    vt = β' * Σx * β + σ2
    mi_true(idx) = (c = Σx[idx, :] * β; 0.5 * log2(vt / (vt - c' * (Σx[idx, idx] \ c))))
    R = min(mi_true([1]), mi_true([2]))
    truth = (R, mi_true([1]) - R, mi_true([2]) - R, mi_true([1, 2]) - mi_true([1]) - mi_true([2]) + R)

    pid = PartialInformationDecomposition(IMMI())
    gauss = two_source_atoms(association(GaussianDecomposition(pid), x1, x2, t))
    @test all(isapprox.(gauss, truth; atol = 0.015))
    # The same closed form from the sample covariance.
    C = cov(hcat(X, t))
    mi_sample(idx) = 0.5 * log2(det(C[idx, idx]) * C[3, 3] / det(C[vcat(idx, 3), vcat(idx, 3)]))
    Rs = min(mi_sample([1]), mi_sample([2]))
    @test gauss[1] ≈ Rs
    @test gauss[4] ≈ mi_sample([1, 2]) - mi_sample([1]) - mi_sample([2]) + Rs
    @test all(isapprox.(two_source_atoms(association(GaussianDecomposition(pid; copula = true),
        exp.(x1), x2 .^ 3, t)), truth; atol = 0.02))
    # `GaussianMI` gives the same closed form, including for the grouped source {12}.
    gmi = two_source_atoms(association(MIDecomposition(pid, GaussianMI()), x1, x2, t))
    @test all(gmi .≈ gauss)

    # Nonparametric estimators, on fewer samples.
    n = 10_000
    for est in (MIDecomposition(pid, KraskovStögbauerGrassberger1(MIShannon(); k = 5)),
            EntropyDecomposition(pid, Kraskov(Shannon(); k = 5)))
        @test all(isapprox.(two_source_atoms(association(est, x1[1:n], x2[1:n], t[1:n])), truth; atol = 0.06))
    end

    # The two discrete routes agree when they bin identically.
    b = CodifyVariables(ValueBinning(RectangularBinning(8)))
    joint = two_source_atoms(association(JointProbabilities(pid, b), x1, x2, t))
    separate = two_source_atoms(association(EntropyDecomposition(pid, PlugIn(Shannon()), b), x1, x2, t))
    @test all(joint .≈ separate)

    # A grouped source: (x1, x2) as one source carries everything, an independent source
    # carries nothing.
    res = association(GaussianDecomposition(pid), StateSpaceSet(X), randn(rng, N), t)
    @test unique_information(res, 1) ≈ mi_true([1, 2]) atol = 0.015
    @test abs(redundancy(res)) < 1e-3 && abs(synergy(res)) < 1e-3
end

# ---------------------------------------------------------------------------------------
# Single atoms, significance testing, accessors and errors.
# ---------------------------------------------------------------------------------------
@testset "PIDAtom and surrogate tests" begin
    N = 2000
    a, b = rand(rng, 0:1, N), rand(rng, 0:1, N)
    t = (a .⊻ b) .⊻ (rand(rng, N) .< 0.1)
    for r in (IMin(), IMMI())
        pid = PartialInformationDecomposition(r)
        res = association(JointProbabilities(pid, disc), a, b, t)
        atom_est(at) = JointProbabilities(PIDAtom(pid, at), disc)
        @test association(atom_est(:synergy), a, b, t) ≈ synergy(res)
        @test association(atom_est(:redundancy), a, b, t) ≈ redundancy(res)
        @test association(atom_est(:unique => 2), a, b, t) ≈ unique_information(res, 2)
        @test association(atom_est(:unique => Int32(2)), a, b, t) ≈ unique_information(res, 2)
        @test association(atom_est([[1], [2]]), a, b, t) ≈ redundancy(res)

        test = SurrogateAssociationTest(atom_est(:synergy); nshuffles = 100, rng = Xoshiro(1))
        @test pvalue(independence(test, a, b, t)) < 0.01
        @test pvalue(independence(test, rand(Xoshiro(5), 0:1, N), b, t)) > 0.01
    end
    # The other estimator types also work through `PIDAtom`.
    x, y = randn(rng, 500), randn(rng, 500)
    z = x .+ y .+ randn(rng, 500)
    pid = PartialInformationDecomposition(IMMI())
    for est in (GaussianDecomposition(pid), MIDecomposition(pid, KraskovStögbauerGrassberger1(MIShannon())),
            EntropyDecomposition(pid, Kraskov(Shannon())),
            EntropyDecomposition(pid, PlugIn(Shannon()), CodifyVariables(ValueBinning(4))))
        atom_est = Associations.with_definition(est, PIDAtom(pid, :synergy))
        @test association(atom_est, x, y, z) ≈ synergy(association(est, x, y, z))
    end
end

@testset "Accessors, printing and errors" begin
    x, y = [0, 0, 1, 1], [0, 1, 0, 1]
    res = association(pid_estimator(IMin()), x, y, x .& y)
    @test first.(atoms(res)) == [[[1], [2]], [[1]], [[2]], [[1, 2]]]
    @test atom(res, ((2,), (1,))) ≈ redundancy(res)
    @test occursin("{12}", sprint(show, MIME("text/plain"), res))

    pid = PartialInformationDecomposition(IMin())
    @test_throws ArgumentError association(GaussianDecomposition(pid), randn(50), randn(50), randn(50))
    @test_throws ArgumentError association(pid_estimator(IMin()), x, x)                  # one source
    @test_throws ArgumentError association(pid_estimator(IMin()), ntuple(_ -> x, 7)...)  # six sources
    @test_throws ArgumentError atom(res, [[1], [1, 2]])  # not an antichain
    @test_throws ArgumentError atom(res, [[3]])          # no such source
    @test_throws ArgumentError PIDAtom(pid, :foo)
    @test_throws ArgumentError PIDAtom(pid, :unique => "a")
end
