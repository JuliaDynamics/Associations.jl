using Statistics: mean, std

function t_statistic(x::AbstractVector; hypothetical_μ = 0.0)
    μ̄ = mean(x)
    σ̄ = std(x)
    degrees_of_freedom = length(x) - 1
    stderr = (σ̄ / degrees_of_freedom)
    t_statistic = (μ̄ - hypothetical_μ) / stderr
end

include("JointDistanceDistributionTest.jl")
include("PredictiveAsymmetryTest.jl")
include("PATest.jl")
