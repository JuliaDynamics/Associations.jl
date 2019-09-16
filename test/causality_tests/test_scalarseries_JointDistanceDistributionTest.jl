x, y = rand(200), rand(200)
jtest = JointDistanceDistributionTest(B = 20)

@test causality(x, y, jtest) isa Vector{T} where T