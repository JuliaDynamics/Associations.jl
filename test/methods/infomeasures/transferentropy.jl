x, y, z = rand(1000), rand(1000), rand(1000)

# Transfer entropy is asymmetric.
pest = SymbolicPermutation(m = 2)
@test transferentropy(TEShannon(), pest, x, y) != transferentropy(TEShannon(), pest, y, x)

est = Lindner( k = 5)
@test transferentropy(est, x, y) isa Real

est = Zhu1(k = 5)
@test_throws DomainError Zhu1(k = 1)
@test transferentropy(est, x, y) isa Real

est = ZhuSingh(k = 5)
@test transferentropy(est, x, y) isa Real

@testset "Convenience" begin
    est = SymbolicTransferEntropy()
    @test transferentropy(est, x, y, z) >= 0.0
    @test transferentropy(TERenyiJizba(), est, x, y, z) isa Real

    est = Hilbert(ValueHistogram(4), target = Amplitude(), source = Phase())
    @test transferentropy(est, x, y, z) >= 0.0
end
