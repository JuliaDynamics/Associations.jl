using Test
using Associations 
using StableRNGs

rng = StableRNG(123)
x, y, z = rand(rng, 30), rand(rng, 30), rand(rng, 30)

X = StateSpaceSet(x)
Y = StateSpaceSet(y)
Z = StateSpaceSet(z)

nshuffles = 2
est_ord = JointProbabilities(CMIShannon(), CodifyVariables(OrdinalPatterns()))
est_vh = JointProbabilities(CMIShannon(), CodifyVariables(ValueHistogram(3)))
est_dp = JointProbabilities(CMIShannon(), CodifyVariables( Dispersion(m = 2)))
lptest_sp = LocalPermutationTest(est_ord; nshuffles, rng)
lptest_vh = LocalPermutationTest(est_vh; nshuffles, rng)
lptest_dp = LocalPermutationTest(est_dp; nshuffles, rng)
@test independence(lptest_sp, x, y, z) isa LocalPermutationTestResult
@test independence(lptest_vh, x, y, z) isa LocalPermutationTestResult
@test independence(lptest_dp, x, y, z) isa LocalPermutationTestResult
@test independence(lptest_sp, X, Y, Z) isa LocalPermutationTestResult
@test independence(lptest_vh, X, Y, Z) isa LocalPermutationTestResult
@test independence(lptest_dp, X, Y, Z) isa LocalPermutationTestResult
