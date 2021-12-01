
x = [0, 0, 1, 0, 0, 1, 0, 0, 1, 0]
y = [0, 0, 0, 1, 0, 0, 1, 0, 0, 1]

@test mean_observed_penchant(x, y, 1, weighted = false) ≈ 1.0
@test mean_observed_penchant(y, x, 1, weighted = false) ≈ 1/7

@test mean_observed_penchant(x, y, 1, weighted = true) ≈ 1.0
@test mean_observed_penchant(y, x, 1, weighted = true) ≈ 3/63