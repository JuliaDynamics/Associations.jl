
import UncertainData: resample 

UT = Union{UncertainValueDataset, UncertainIndexDataset, UncertainDataset}

function s_measure(source::UT, target::UT, m::Int, τ, K::Int; kwargs...)
    s_measure(resample(source), resample(target), m, τ, K; kwargs...)
end

function jdd(source::UT, target::UT, m::Int, τ, K::Int; kwargs...)
    jdd(resample(source), resample(target), m, τ, K; kwargs...)
end

function transferentropy(source::UT, target::UT, emb, estimator)
    transferentropy(resample(source), resample(target), emb, estimator)
end

function transferentropy(source::UT, target::UT, cond::UT, emb, estimator)
    transferentropy(resample(source), resample(target), resample(cond), emb, estimator)
end

function predictive_asymmetry(source::UT, target::UT, ηs, estimator; kwargs...)
    predictive_asymmetry(resample(source), resample(target), ηs, estimator; kwargs...)
end

function predictive_asymmetry(source::UT, target::UT, cond::UT, ηs, estimator; kwargs...)
    predictive_asymmetry(resample(source), resample(target), resample(cond), ηs, estimator; kwargs...)
end

