
import UncertainData: resample 


##################################################
# Basic resampling for `UncertainDataset`s
##################################################
UT = Union{UncertainValueDataset, UncertainIndexDataset, UncertainDataset}

s_measure(source::UT, target::UT, m::Int, τ, K::Int; kwargs...) =
    s_measure(resample(source), resample(target), m, τ, K; kwargs...)

jdd(source::UT, target::UT, m::Int, τ, K::Int; kwargs...) =
    jdd(resample(source), resample(target), m, τ, K; kwargs...)

mutualinfo(source::UT, target::UT, estimator; kwargs...) =
    mutualinfo(resample(source), resample(target), estimator; kwargs...)

transferentropy(source::UT, target::UT, estimator; kwargs...) =
    transferentropy(resample(source), resample(target), estimator; kwargs...)

transferentropy(source::UT, target::UT, cond::UT, estimator; kwargs...) =
    transferentropy(resample(source), resample(target), resample(cond), estimator; kwargs...)

predictive_asymmetry(source::UT, target::UT, estimator; kwargs...) =
    predictive_asymmetry(resample(source), resample(target), estimator; kwargs...)

predictive_asymmetry(source::UT, target::UT, cond::UT, estimator; kwargs...) =
    predictive_asymmetry(resample(source), resample(target), resample(cond), estimator; kwargs...)

crossmap(source::UT, target::UT; kwargs...) = 
    crossmap(resample(source), resample(target); kwargs...)

ccm(source::UT, target::UT, timeseries_lengths; kwargs...) = 
    ccm(resample(source), resample(target), timeseries_lengths; kwargs...)

##########################################################################
# Basic resampling for `UncertainIndexValueDataset` (no constraints)
##########################################################################
UIVD = UncertainIndexValueDataset

# TODO: warn about potential index reversals?
#
# function warn_about_sampling(source::V, target::W) 
#     if source isa UIVD
#         @warn "`source` isa UncertainIndexValueDataset. Index reversals may occur. Consider constrained resampling."
#     end

#     if target isa UIVD
#         @warn "`target` isa UncertainIndexValueDataset. Index reversals may occur. Consider constrained resampling."
#     end
# end

# function warn_about_sampling(source::V, target::W, cond::X)
#     warn_about_sampling(source, target)
#     if cond isa UIVD
#         @warn "`cond` isa UncertainIndexValueDataset. Index reversals may occur. Consider constrained resampling."
#     end
# end

s_measure(source::UIVD, target::UIVD, m::Int, τ, K::Int; kwargs...)
    s_measure(resample(source.values), resample(target.values), m, τ, K; kwargs...)

jdd(source::UIVD, target::UIVD, m::Int, τ, K::Int; kwargs...) =
    jdd(resample(source.values), resample(target.values), m, τ, K; kwargs...)

mutualinfo(source::UIVD, target::UIVD, estimator; kwargs...) =
    mutualinfo(resample(source.values), resample(target.values), estimator; kwargs...)

transferentropy(source::UIVD, target::UIVD, estimator; kwargs...) =
    transferentropy(resample(source.values), resample(target.values), estimator; kwargs...)

transferentropy(source::UIVD, target::UIVD, cond::UIVD, estimator; kwargs...) =
    transferentropy(resample(source.values), resample(target.values), resample(cond.values), estimator; kwargs...)

predictive_asymmetry(source::UIVD, target::UIVD, estimator; kwargs...) =
    predictive_asymmetry(resample(source.values), resample(target.values), estimator; kwargs...)

predictive_asymmetry(source::UIVD, target::UIVD, cond::UIVD, estimator; kwargs...) =
    predictive_asymmetry(resample(source.values), resample(target.values), resample(cond.values), estimator; kwargs...)

crossmap(source::UIVD, target::UIVD; kwargs...) = 
    crossmap(resample(source.values), resample(target.values); kwargs...)

ccm(source::UIVD, target::UIVD, timeseries_lengths; kwargs...) = 
    ccm(resample(source.values), resample(target.values), timeseries_lengths; kwargs...)

