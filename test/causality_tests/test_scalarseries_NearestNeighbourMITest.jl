import StatsBase

x, y = rand(100), rand(100)
@test causality(x, y, NearestNeighbourMITest(ηs = 5)) isa Vector
@test causality(x, y, NearestNeighbourMITest(ηs = 2:10)) isa Vector
