
import UncertainData: resample 

UT = Union{UncertainValueDataset, UncertainIndexDataset, UncertainDataset}

function s_measure(source::UT, target::UT, m::Int, τ, K::Int; kwargs...)
    s_measure(resample(source), resample(target), m, τ, K; kwargs...)
end

function jdd(source::UT, target::UT, m::Int, τ, K::Int; kwargs...)
    jdd(resample(source), resample(target), m, τ, K; kwargs...)
end

function transferentropy(source::UT, target::UT, estimator; kwargs...)
    transferentropy(resample(source), resample(target), estimator; kwargs...)
end

function transferentropy(source::UT, target::UT, cond::UT; kwargs...)
    transferentropy(resample(source), resample(target), resample(cond), estimator; kwargs...)
end

function predictive_asymmetry(source::UT, target::UT, estimator; kwargs...)
    predictive_asymmetry(resample(source), resample(target), estimator; kwargs...)
end

function predictive_asymmetry(source::UT, target::UT, cond::UT, estimator; kwargs...)
    predictive_asymmetry(resample(source), resample(target), resample(cond), estimator; kwargs...)
end

