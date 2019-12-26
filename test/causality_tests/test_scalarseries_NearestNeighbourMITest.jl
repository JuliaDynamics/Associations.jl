import StatsBase

x, y = rand(100), rand(100)

####################
# Single binnings
####################

@test causality(x, y, NearestNeighbourMITest(ηs = 5)) isa Array{Float64,0}
@test causality(x, y, NearestNeighbourMITest(ηs = 2:10)) isa Vector
