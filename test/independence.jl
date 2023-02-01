x, y, z = rand(100), rand(100), rand(100)

@testset "Conditional independence tests" begin
    test_cmi = LocalPermutation(CMIShannon(), FPVP())
    @test independence(test_cmi, x, y, z) isa LocalPermutationTestResult

    # Not yet implemented
    test_te = LocalPermutation(TEShannon(), FPVP())
    @test_throws ArgumentError independence(test_te, x, y, z) isa LocalPermutationTestResult

end
