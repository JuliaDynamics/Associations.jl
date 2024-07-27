using Test
using Associations 
using StableRNGs

rng = StableRNG(123)
x, y, z = rand(rng, 30), rand(rng, 30), rand(rng, 30)

independence_test = LocalPermutationTest(FPVP(CMIShannon()))
# We should get back a convenience wrapper containing the result.
res = independence(independence_test, x, z, y)
@test res isa LocalPermutationTestResult

# We should be able to compute p-values for the result.
@test pvalue(res) isa Real
@test pvalue(res) ≥ 0

# Only conditional analyses are possible, meaning that we need three inputs.
# Pairwise analyses won't work, because only two inputs are given.
@test_throws ArgumentError independence(independence_test, x, y)

# Sampling with/without replacement
test_cmi_replace = LocalPermutationTest(FPVP(CMIShannon()), replace = true)
test_cmi_nonreplace = LocalPermutationTest(FPVP(CMIShannon()), replace = false)
@test independence(test_cmi_replace, x, y, z) isa LocalPermutationTestResult
@test independence(test_cmi_nonreplace, x, y, z) isa LocalPermutationTestResult

# Measure definition AND estimator must be provided for info measures
@test_throws ArgumentError LocalPermutationTest(TEShannon()) # estimator needed

# The number of local neighbors can't exceed the number of input datapoints
test_kperm_toolarge = LocalPermutationTest(FPVP(CMIShannon()); kperm = 200, rng)
@test_throws ArgumentError independence(test_kperm_toolarge, x, y, z)


# --------------------------------
# Output
# --------------------------------
rng = StableRNG(123)
x, y, z = rand(rng, 30), rand(rng, 30), rand(rng, 30)
independence_test = LocalPermutationTest(FPVP(CMIShannon()), nshuffles = 2)
res = independence(independence_test, x, z, y)

# Internals
out_str_pval = repr( Associations.pvalue_text_summary(res))
@test occursin("p-value:", out_str_pval)

# Internals
out_str_pval = repr( Associations.quantiles_text(res))
@test occursin("Ensemble quantiles", out_str_pval)


out_str_conclusion = repr( Associations.null_hypothesis_text(res))
@test occursin("The first two variables are independent, given the 3rd variable", out_str_conclusion)

independence_test = SurrogateAssociationTest(KSG1(MIShannon()))
res = independence(independence_test, x, y)
out_str_conclusion = repr( Associations.null_hypothesis_text(res))
@test occursin("The variables are independent", out_str_conclusion)
