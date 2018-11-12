d1 = rand(50)
d2 = rand(50)
r = rand(50)

# Directly from time series
@test te(d1, r, estimator = :tetogrid) >= 0
@test te(d2, r, estimator = :tetogrid) >= 0
@test te(d1, r, estimator = :tefreq) >= 0
@test te(d2, r, estimator = :tefreq) >= 0
@test te(d1, r, estimator = :tekNN) >= 0
@test te(d2, r, estimator = :tekraskov) >= 0

#########################
# Surrogates
#########################
@test te(d1, r, which_is_surr = 1) >= 0
@test te(d1, r, which_is_surr = 2) >= 0

@test te(d1, r, which_is_surr = 2, surr_type = aaft) >= 0
@test te(d1, r, which_is_surr = 2, surr_type = iaaft) >= 0
@test te(d1, r, which_is_surr = 2, surr_type = randomphases) >= 0
@test te(d1, r, which_is_surr = 2, surr_type = randomamplitudes) >= 0
@test te(d1, r, which_is_surr = 2, surr_type = randomshuffle) >= 0

#########################
# Forward prediction lags
#########################
@test te(d1, r, ν = 1) >= 0
@test te(d1, r, ν = 5) >= 0 

#########################
# Embedding lags and dimensions
#########################
d1 = rand(100)
d2 = rand(100)
r = rand(100)
@test te(d1, r, ν = 1, τ = 2, dim = 4) >= 0
@test te(d1, r, ν = 3, τ = 3, dim = 5) >= 0

#########################
# Tuning the bins
#########################
@test te(d1, r, n_ϵ = 10) >= 0
@test te(d1, r, n_ϵ = 10, max_numbins = 5, min_numbins = 2) >= 0



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
