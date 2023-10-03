using StableRNGs
rng = StableRNG(123)
x, y, z = rand(rng, 100), rand(rng, 100), rand(rng, 100)

test_cmi_replace = LocalPermutationTest(CMIShannon(), FPVP())
test_cmi_nonreplace = LocalPermutationTest(CMIShannon(), FPVP())

test_teshannon = LocalPermutationTest(TEShannon(), FPVP())
@test_throws ArgumentError LocalPermutationTest(TEShannon()) # estimator needed

@test independence(test_cmi_replace, x, y, z) isa LocalPermutationTestResult
@test independence(test_cmi_nonreplace, x, y, z) isa LocalPermutationTestResult

@test independence(test_teshannon, x, y, z) isa LocalPermutationTestResult

test_kperm_toolarge = LocalPermutationTest(CMIShannon(), FPVP(); kperm = 200, rng)
@test_throws ArgumentError independence(test_kperm_toolarge, x, y, z)

# CMI
# ------------------------
# Independence tests
x = rand(rng, 50)
y = rand(rng, 50)
z = rand(rng, 50)
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

# Part mutual information
# ------------------------
# Independence tests
x = rand(rng, 50)
y = rand(rng, 50)
z = rand(rng, 50)
X = StateSpaceSet(x)
Y = StateSpaceSet(y)
Z = StateSpaceSet(z)

nshuffles = 5
lptest_sp = LocalPermutationTest(PMI(), SymbolicPermutation(); nshuffles, rng)
lptest_vh = LocalPermutationTest(PMI(), ValueHistogram(4); nshuffles, rng)
lptest_dp = LocalPermutationTest(PMI(), Dispersion(); nshuffles, rng)
@test independence(lptest_sp, x, y, z) isa LocalPermutationTestResult
@test independence(lptest_vh, x, y, z) isa LocalPermutationTestResult
@test independence(lptest_dp, x, y, z) isa LocalPermutationTestResult
@test independence(lptest_sp, X, Y, Z) isa LocalPermutationTestResult
@test independence(lptest_vh, X, Y, Z) isa LocalPermutationTestResult
@test independence(lptest_dp, X, Y, Z) isa LocalPermutationTestResult


α = 0.05
n = 10000
x = rand(rng, 1:3, n);
y = map(x) do xᵢ
    if xᵢ == 1
        return rand(rng, 1:2)
    elseif xᵢ == 21
        return rand(rng, 2:3)
    else
        return rand(rng, 4:5)
    end
end;
z = map(y) do yᵢ
    if yᵢ == 1 || yᵢ == 2
        return rand(rng, [1, 3])
    elseif yᵢ == 2 || yᵢ == 3
        return rand(rng, [2, 4])
    else
        return rand(rng, [5, 7])
    end
end;
# Add some noise, so that local permutation test works
x = x + rand(rng, n) * 1e-3
y = y + rand(rng, n) * 1e-3
z = z + rand(rng, n) * 1e-3

# We should not be able to reject the null hypothesis `x ⫫ z | y`, because
# x → y → z, so when conditioning on the intermediate variable,
# the first and last variable in the chain should be independent.
test_sp = LocalPermutationTest(PMI(), SymbolicPermutation(); nshuffles = 200, rng)
test_dp = LocalPermutationTest(PMI(), Dispersion(); nshuffles = 200, rng)
test_vh = LocalPermutationTest(PMI(), ValueHistogram(2); nshuffles = 200, rng)
@test pvalue(independence(test_sp, x, y, z)) > α
@test pvalue(independence(test_dp, x, y, z)) > α
@test pvalue(independence(test_vh, x, y, z)) > α

# tests for specific measures
include("LocalPermutationTest/transferentropy.jl")
