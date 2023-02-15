x, y, z = rand(100), rand(100), rand(100)

test_cmi_replace = LocalPermutationTest(CMIShannon(), FPVP())
test_cmi_nonreplace = LocalPermutationTest(CMIShannon(), FPVP())

test_teshannon = LocalPermutationTest(TEShannon(), FPVP())
@test_throws ArgumentError LocalPermutationTest(TEShannon()) # estimator needed

@test independence(test_cmi_replace, x, y, z) isa LocalPermutationTestResult
@test independence(test_cmi_nonreplace, x, y, z) isa LocalPermutationTestResult

@test independence(test_teshannon, x, y, z) isa LocalPermutationTestResult

test_kperm_toolarge = LocalPermutationTest(CMIShannon(), FPVP(), kperm = 200)
@test_throws ArgumentError independence(test_kperm_toolarge, x, y, z)

# # Not yet implemented
# test_te = LocalPermutationTest(TEShannon(), FPVP())
# @test_throws ArgumentError independence(test_te, x, y, z) isa LocalPermutationTestResult
