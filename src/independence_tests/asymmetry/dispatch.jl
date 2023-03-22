function dispatch(measure, est, x, y)
    estimate(measure, est, x, y)
end

function dispatch(measure, est, x, y, z)
    estimate(measure, est, x, y, z)
end

const CMI_MEASURES = Union{
    ConditionalMutualInformation,
    MutualInformation}

function dispatch(measure::CMI_MEASURES, est, x, y, z)
    condmutualinfo(measure, est, x, y, z)
end

function dispatch(measure::CMIShannon, est, x, y)
    mutualinfo(MIShannon(; base = measure.base), est, x, y)
end
