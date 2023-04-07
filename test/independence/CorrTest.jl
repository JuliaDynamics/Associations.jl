using HypothesisTests: CorrelationTest, pvalue
using StableRNGs
rng = StableRNG(123)
X = randn(rng, 500);
Y = randn(rng, 500);
Z = randn(rng, 500) .+ Y;

pr = independence(CorrTest(), X, Y);
prc = independence(CorrTest(), X, Y);
@test -1.0 ≤ independence(CorrTest(), X, X).ρ ≤ 1.0
@test -1.0 ≤ independence(CorrTest(), X, -X).ρ ≤ 1.0
@test -1.0 ≤ independence(CorrTest(), X, X, Y).ρ ≤ 1.0
@test -1.0 ≤ independence(CorrTest(), X, -X, Y).ρ ≤ 1.0
@test pvalue(pr) ≥ 0.0
@test pvalue(prc) ≥ 0.0
# ----------------------------------------------------------------
# "Analytical" tests: compare to HypothesisTests.jl.
# ----------------------------------------------------------------
prh = CorrelationTest(X, Y);
prhc = CorrelationTest(X, Y);

# Computed coefficients are roughly equal.
@test prh.r ≈ pr.ρ
@test prhc.r ≈ prc.ρ

# Computed p-value are roughly equal
@test round(pvalue(pr), digits = 3) ≈ round(pvalue(prh), digits = 3)
@test round(pvalue(prc), digits = 3) ≈ round(pvalue(prhc), digits = 3)
