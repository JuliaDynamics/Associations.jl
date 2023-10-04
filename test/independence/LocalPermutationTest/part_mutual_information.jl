using Test
using CausalityTools 
using StableRNGs

rng = StableRNG(123)
x, y, z = rand(rng, 30), rand(rng, 30), rand(rng, 30)

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
