d1 = rand(50)
d2 = rand(50)
r = rand(50)

# Directly from time series
@test te(d1, r) >= 0
@test te(d2, r) >= 0
# maxlag = 2
# res1 = te_lags(d1, r, -maxlag:maxlag, method = :transferoperator_grid,
#     max_numbins = 7, min_numbins = 3)
# res2 = te_lags(d1, r, -maxlag:maxlag, method = :visitfreq,
#     max_numbins = 10, min_numbins = 3)
# res3 = te_lags(d1, r, -maxlag:maxlag, method = :transferoperator_grid)
# res4 = te_lags(d1, r, -maxlag:maxlag, method = :visitfreq)
# res5 = te_lags(d1, r, -maxlag:maxlag, method = :kraskov)
#
# @test all(res1 .>= 0)
# @test all(res2 .>= 0)
# @test all(res3 .>= 0)
# @test all(res4 .>= 0)
# @test all(res5 .>= 0)
#
# # From embedding with associated TEVars instance
# E = embed([d1, d2, r], [1, 2, 3, 3], [0, 0, 0, 1])
# v = TEVars([4], [3], [1, 2])
# res1 = te_lags(E, v, -maxlag:maxlag, method = :transferoperator_grid,
#     max_numbins = 7, min_numbins = 3)
# res2 = te_lags(E, v, -maxlag:maxlag, method = :visitfreq,
#     max_numbins = 10, min_numbins = 3)
# res3 = te_lags(E, v, -maxlag:maxlag, method = :transferoperator_grid)
# res4 = te_lags(E, v, -maxlag:maxlag, method = :visitfreq)
# res5 = te_lags(E, v, -maxlag:maxlag, method = :kraskov)
#
# @test all(res1 .>= 0)
# @test all(res2 .>= 0)
# @test all(res3 .>= 0)
# @test all(res4 .>= 0)
# @test all(res5 .>= 0)
