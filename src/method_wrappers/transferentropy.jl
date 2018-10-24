include("TE.jl")

# """
#     te_lags(driver, response, lags::Union{Vector{Int}, UnitRange{Int}};
#            method = :transferoperator_grid,
#            E = nothing, v = nothing,
#            ϵ = nothing, n_ϵ = 5, dim = 3, η = 1,
#            k1 = 4, k2 = 5, metric = Chebyshev(),
#            which_is_surr = -1, surr_type = aaft)
# """
# function te_lags(driver, response, lags::Union{Vector{Int}, UnitRange{Int}};
#         method = :transferoperator_grid,
#         E = nothing, v = nothing,
#         ϵ = nothing, n_ϵ = 5, dim = 3, η = 1,
#         k1 = 4, k2 = 5, metric = Chebyshev(),
#         which_is_surr = -1, surr_type = aaft)
#
#     pmap(l -> te(driver, response,
#                     lag = l,
#                     method = method,
#                     E = E, v = v,
#                     ϵ = ϵ, n_ϵ = n_ϵ, dim = dim, η = η,
#                     k1 = k1, k2 = k2, metric = metric,
#                     which_is_surr = which_is_surr,
#                     surr_type = surr_type), lags)
# end
#
# """
#     te_surr(driver, response, nsurr::Int; lag = 1,
#             method = :transferoperator_grid,
#             E = nothing, v = nothing,
#             ϵ = nothing, n_ϵ = 5, dim = 3, η = 1,
#             k1 = 4, k2 = 5, metric = Chebyshev(),
#             which_is_surr = 0, surr_type = aaft)
# """
# function te_surr(driver, response, nsurr::Int; lag = 1,
#         method = :transferoperator_grid,
#         E = nothing, v = nothing,
#         ϵ = nothing, n_ϵ = 5, dim = 3, η = 1,
#         k1 = 4, k2 = 5, metric = Chebyshev(),
#         which_is_surr = 0, surr_type = aaft)
#
#     pmap(l -> te(driver, response,
#                     lag = lag,
#                     method = method,
#                     E = E, v = v,
#                     ϵ = ϵ, n_ϵ = n_ϵ, dim = dim, η = η,
#                     k1 = k1, k2 = k2, metric = metric,
#                     which_is_surr = which_is_surr,
#                     surr_type = surr_type), 1:nsurr)
# end

export te#, te_lags, te_surr

#
# function te_surrnormalized_p(driver, response; ϵ = 6, lag = 1, surr_type = aaft,
#                             which_is_surr = 1, n_surr = 10, q = 0.95)
#     te = te_original(driver, response, ϵ = ϵ, lag = lag)
#     surr = SharedVector{Float64}(n_surr)
#     @sync @parallel for i = 1:n_surr
#         surr[i] = te_surr(driver, response, surr_type = surr_type,
#             ϵ = ϵ, lag = lag, which_is_surr = which_is_surr)
#     end
#     te - quantile(Array(surr), q)
# end
#
#
# function te_surrnormalized(driver, response; ϵ = 6, lag = 1, surr_type = aaft,
#                             which_is_surr = 1, n_surr = 10, q = 0.95)
#     te = te_original(driver, response, ϵ = ϵ, lag = lag)
#     surr = [te_surr(driver, response, surr_type = surr_type,
#             ϵ = ϵ, lag = lag, which_is_surr = which_is_surr) for i = 1:n_surr]
#     te - quantile(surr, q)
# end
#
# function te_regular_and_surr(driver, response; ϵ = 6, lag = 1, surr_type = aaft,
#                             which_is_surr = 1, n_surr = 10, q = 0.95)
#     te = te_original(driver, response, ϵ = ϵ, lag = lag)
#     surr = [te_surr(driver, response, surr_type = surr_type,
#             ϵ = ϵ, lag = lag, which_is_surr = which_is_surr) for i = 1:n_surr]
#     te, surr
# end
#

#
#
# function te_auto_lags_wsurr(driver, response, lags)
#     pmap(l -> te_binsizes(driver, response, lag = l), lags)
# end
#
#
# export te_original,
#     te_original_lags
#     te_binsizes,
#     te_surr,
#     te_surrnormalized,
#     te_surrnormalized_p,
#     te_regular_and_surr,
#     te_auto_lags,
#     te_auto_lags_wsurr
