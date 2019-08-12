# using Test, Distances, TransferEntropy

# d1 = rand(50)
# d2 = rand(50)
# r = rand(50)
# tol = 1e-8

# # Directly from time series
# @test standard_te(d1, r, estimator = :tetogrid) >= 0 - tol
# @test standard_te(d2, r, estimator = :tetogrid) >= 0 - tol
# @test standard_te(d1, r, estimator = :tefreq) >= 0 - tol
# @test standard_te(d2, r, estimator = :tefreq) >= 0 - tol
# @test standard_te(d1, r, estimator = :tekNN) >= 0 - tol
# @test standard_te(d2, r, estimator = :tekraskov) >= 0 - tol

# #########################
# # Surrogates
# #########################
# @test standard_te(d1, r, which_is_surr = :driver) >= 0 - tol
# @test standard_te(d1, r, which_is_surr = :response) >= 0 - tol
# @test standard_te(d1, r, which_is_surr = :none) >= 0 - tol

# @test standard_te(d1, r, which_is_surr = :both, surr_func = aaft) >= 0 - tol
# @test standard_te(d1, r, which_is_surr = :both, surr_func = iaaft) >= 0 - tol
# @test standard_te(d1, r, which_is_surr = :both, surr_func = randomphases) >= 0 - tol
# @test standard_te(d1, r, which_is_surr = :both, surr_func = randomamplitudes) >= 0 - tol
# @test standard_te(d1, r, which_is_surr = :both, surr_func = randomshuffle) >= 0 - tol

# #########################
# # Forward prediction lags
# #########################
# @test standard_te(d1, r, ν = 1) >= 0 - tol
# @test standard_te(d1, r, ν = 5) >= 0 - tol

# #########################
# # Embedding lags and dimensions
# #########################
# d1 = rand(100)
# d2 = rand(100)
# r = rand(100)
# @test standard_te(d1, r, ν = 1, τ = 2, dim = 4) >= 0 - tol
# @test standard_te(d1, r, ν = 3, τ = 3, dim = 5) >= 0 - tol

# #########################
# # Tuning the bins
# #########################
# @test standard_te(d1, r, n_ϵ = 10) >= 0 - tol
# @test standard_te(d1, r, n_ϵ = 10, max_numbins = 5, min_numbins = 2) >= 0 - tol

##########################
# Regular TE
##########################

using Test, TransferEntropy, CausalityToolsBase
k, l, m = 1, 1, 1
τ = 1
η = 1
n_subdivs = 3
te_res1 = te_reg(rand(1000), rand(1000), k, l, m, η = η, τ = τ, n_subdivs = n_subdivs) 
te_res2 = te_reg(rand(1000), rand(1000), k + 1, l, m + 1, η = η, τ = τ + 1, n_subdivs = n_subdivs) 

@test te_res1 isa Vector{<:Real}
@test te_res2 isa Vector{<:Real}
@test length(te_res1) == 4


@test te_reg(rand(100), rand(100), 1, 1, 1, RectangularBinning(4)) |> typeof <: Real
@test te_reg(rand(100), rand(100), 1, 1, 1, [RectangularBinning(4), RectangularBinning(0.2)]) isa Vector{<:Real}

##########################
# Conditional TE
##########################

# Default discretization schemes
using Test
k, l, m, n = 1, 1, 1, 1
τ = 1
η = 1
n_subdivs = 3
te_res1 = te_cond(rand(1000), rand(1000), rand(1000), k, l, m, n, η = η, τ = τ, n_subdivs = n_subdivs) 
te_res2 = te_cond(rand(1000), rand(1000), rand(1000), k + 1, l, m + 1, n, η = η, τ = τ + 1, n_subdivs = n_subdivs) 

@test te_res1 isa Vector{<:Real}
@test te_res2 isa Vector{<:Real}
@test length(te_res1) == 4


# User-provided binning schemes
using Test
x, y, z = rand(100), rand(100), rand(100)
@test te_cond(x, y, z, 1, 1, 1, 1, RectangularBinning(4)) |> typeof <: Real
@test te_cond(x, y, z, 1, 1, 1, 1, [RectangularBinning(4), RectangularBinning(0.2)]) isa Vector{<:Real}
@test te_cond(x, y, z, 1, 1, 1, 1, RectangularBinning(3), estimator = TransferOperatorGrid()) |> typeof <: Real