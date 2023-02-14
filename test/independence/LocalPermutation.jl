x, y, z = rand(100), rand(100), rand(100)

test_cmi = LocalPermutationTest(CMIShannon(), FPVP())
@test independence(test_cmi, x, y, z) isa LocalPermutationTestResult

# # Not yet implemented
# test_te = LocalPermutationTest(TEShannon(), FPVP())
# @test_throws ArgumentError independence(test_te, x, y, z) isa LocalPermutationTestResult
