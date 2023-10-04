using Test
using CausalityTools 
using StableRNGs

rng = StableRNG(123)
x, y, z = rand(rng, 30), rand(rng, 30), rand(rng, 30)

X = StateSpaceSet(x)
Y = StateSpaceSet(y)
Z = StateSpaceSet(z)

nshuffles = 5
lptest_sp = LocalPermutationTest(CMIShannon(), SymbolicPermutation(); nshuffles, rng)
lptest_vh = LocalPermutationTest(CMIShannon(), ValueHistogram(4); nshuffles, rng)
lptest_dp = LocalPermutationTest(CMIShannon(), Dispersion(); nshuffles, rng)
@test independence(lptest_sp, x, y, z) isa LocalPermutationTestResult
@test independence(lptest_vh, x, y, z) isa LocalPermutationTestResult
@test independence(lptest_dp, x, y, z) isa LocalPermutationTestResult
@test independence(lptest_sp, X, Y, Z) isa LocalPermutationTestResult
@test independence(lptest_vh, X, Y, Z) isa LocalPermutationTestResult
@test independence(lptest_dp, X, Y, Z) isa LocalPermutationTestResult
