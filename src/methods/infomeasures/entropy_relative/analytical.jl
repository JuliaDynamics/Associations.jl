# This function contain analytical expressions for various relative entropies.
# Renyi divergence expressions are from Gil, M. (2011). On Rényi divergence measures for continuous alphabet sources. PhD Thesis.
using Distributions: Beta
using SpecialFunctions: gamma
_beta(x, y) = gamma(x)*gamma(y) / gamma(x + y)

function entropy_relative(e::Renyi, x::Beta, y::Beta)
    q = e.q # Our q is their α
    αx, αy, βx, βy = x.α, y.α, x.β, y.β
    a = q*αx + (1 - q)*αy
    b = q*βx + (1 - q)*βy
    @assert a >= 0 && b >= 0
    re = log(_beta(αy, βy) / _beta(αx, βx)) +
        (1 / (q - 1)) * log(_beta(a, b) / _beta(αx, βx))
    return re / log(e.base, ℯ)
end
