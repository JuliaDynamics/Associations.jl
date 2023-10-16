using Statistics: mean, std

function t_statistic(x::AbstractVector; hypothetical_μ = 0.0)
    μ̄ = mean(x)
    σ̄ = std(x)
    degrees_of_freedom = length(x) - 1
    stderr = (σ̄ / degrees_of_freedom)
    t_statistic = (μ̄ - hypothetical_μ) / stderr
end


"""
    fishers_z(p̂)

Compute Fisher's z-transform on the sample partial correlation coefficient `p̂` (computed
as the correlation between variables `i` and `j`, given the remaining variables):

```math
Z(V_i, V_j | \\bf{S}) = \\dfrac{1}{2}
\\log{\\left(
    \\dfrac{1 + \\hat{p}(V_i, V_j | \\bf{S})}{1 - \\hat{p}(V_i, V_j | \\bf{S})}
\\right) }
```
"""
function fishers_z(p̂)
    return 0.5 * log((1 + p̂) / (1 - p̂))
end

include("JointDistanceDistributionTest.jl")
#include("PredictiveAsymmetryTest.jl")
#include("PATest.jl")
#include("CorrTest.jl")
