using Test, Distances
@testset "TE wrapper" begin
    d1 = rand(50)
    d2 = rand(50)
    r = rand(50)
    tol = 1e-8

    # Directly from time series
    te(d1, r, estimator = :tetogrid)
    @test te(d1, r, estimator = :tetogrid) >= 0 - tol
    @test te(d2, r, estimator = :tetogrid) >= 0 - tol
    @test te(d1, r, estimator = :tefreq) >= 0 - tol
    @test te(d2, r, estimator = :tefreq) >= 0 - tol
    @test te(d1, r, estimator = :tekNN) >= 0 - tol
    @test te(d2, r, estimator = :tekraskov) >= 0 - tol

    #########################
    # Surrogates
    #########################
    @test te(d1, r, which_is_surr = :driver) >= 0 - tol
    @test te(d1, r, which_is_surr = :response) >= 0 - tol
    @test te(d1, r, which_is_surr = :none) >= 0 - tol

    @test te(d1, r, which_is_surr = :both, surr_func = aaft) >= 0 - tol
    @test te(d1, r, which_is_surr = :both, surr_func = iaaft) >= 0 - tol
    @test te(d1, r, which_is_surr = :both, surr_func = randomphases) >= 0 - tol
    @test te(d1, r, which_is_surr = :both, surr_func = randomamplitudes) >= 0 - tol
    @test te(d1, r, which_is_surr = :both, surr_func = randomshuffle) >= 0 - tol

    #########################
    # Forward prediction lags
    #########################
    @test te(d1, r, ν = 1) >= 0 - tol
    @test te(d1, r, ν = 5) >= 0 - tol

    #########################
    # Embedding lags and dimensions
    #########################
    d1 = rand(100)
    d2 = rand(100)
    r = rand(100)
    @test te(d1, r, ν = 1, τ = 2, dim = 4) >= 0 - tol
    @test te(d1, r, ν = 3, τ = 3, dim = 5) >= 0 - tol

    #########################
    # Tuning the bins
    #########################
    @test te(d1, r, n_ϵ = 10) >= 0 - tol
    @test te(d1, r, n_ϵ = 10, max_numbins = 5, min_numbins = 2) >= 0 - tol

end
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
