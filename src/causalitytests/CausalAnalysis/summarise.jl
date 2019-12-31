const PT{T, N} = Union{PredictiveAsymmetryTest{T, N}, 
    NormalisedPredictiveAsymmetryTest{T, N}} where {T, N}
const SPT{T, N} = SubsampledDataCausalityTest{<:PT{T, N}} where {T, N}

"""
    summarise(f::Function, ::CausalAnalysis{MetaCausalityTest}, args...; kwargs...)

Summarise a causality analysis performed by repeated subsampling of the input data using 
the statistic `f`. The statistic is apply on a per-subsample basis, so that one value 
for the summary statistic is returned for each evaluation of the causality test.
"""
function summarise(f::Function, analysis::CausalAnalysis{MetaCausalityTest}, 
        args...; kwargs...)
    f.(analysis.result, args...; kwargs...)
end

""" 
For predictive tests, don't summarise on a per-subsample basis, but gather the 
subsamples and summarise per prediction lag (i.e. over the `ηs`). 
"""
function summarise(f::Function, analysis::CausalAnalysis{CT, DT, RT}, 
        args...; kwargs...) where {CT <: SPT{T, N}, DT, RT} where {T, N}

    M = hcat(analysis.result...)
    η_results = [M[i, :] for i = 1:N]
    f.(η_results, args...; kwargs...)
end