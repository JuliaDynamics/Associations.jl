x, y = rand(100), rand(100)
X, Y = Dataset(rand(100, 3)), Dataset(rand(100, 2))
Z, W = Dataset(rand(110, 2)), Dataset(rand(90, 4))
dx, τx = 2, 1
dy, τy = 2, 1

# V2.X
@test h_measure(HMeasure(), x, y) isa Float64
@test h_measure(HMeasure(dx = dx, τx = τx), x, Y) isa Float64
@test h_measure(HMeasure(dy = dy, τy = τy), X, y) isa Float64
@test h_measure(HMeasure(), X, Y) isa Float64
# test that multivariate datasets are being length-matched
@test h_measure(HMeasure(), X, Z) isa Float64
@test h_measure(HMeasure(), W, X) isa Float64
