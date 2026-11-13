using Test
using Associations
using Random
rng = Xoshiro(1234)

# `codified_marginals` with `ValueBinning(RectangularBinning(…))` gives one encoding per
# input, also when an input is a multi-column `StateSpaceSet`.
@testset "Multi-column inputs with ValueBinning" begin
    N = 20_000
    a, b, c = randn(rng, N), randn(rng, N), randn(rng, N)
    y = a .+ 0.1 .* randn(rng, N)
    G = StateSpaceSet(a, b) # a two-column input
    d = CodifyVariables(ValueBinning(RectangularBinning(4)))

    enc = codified_marginals(d, G, y)
    @test length(enc) == 2
    @test dimension(enc[1]) == 2
    @test enc[2] isa AbstractVector
    # Each encoding holds the bins of its own input's columns.
    enc_a, enc_b, enc_y = codified_marginals(d, a, b, y)
    @test Matrix(enc[1]) == hcat(enc_a, enc_b)
    @test enc[2] == enc_y

    # y is nearly a copy of the first column of G, so I(G; y) is large. (Before the fix,
    # y was replaced by the second column of G, which is independent of the first.)
    mi_sep = association(EntropyDecomposition(MIShannon(), PlugIn(Shannon()), d), G, y)
    mi_joint = association(JointProbabilities(MIShannon(), d), G, y)
    @test mi_sep > 0.5
    @test mi_sep ≈ mi_joint

    # c is independent of everything, so I(G; y | c) ≈ I(G; y), and I(G; c | y) ≈ 0.
    cmi_est = EntropyDecomposition(CMIShannon(), PlugIn(Shannon()), d)
    @test association(cmi_est, G, y, c) > 0.5
    @test association(cmi_est, G, c, y) < 0.05
end
