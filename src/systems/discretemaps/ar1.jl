"""
    eom_ar1(x, p, n) -> Function

Equations of motion for a bivariate, order one autoregressive model, where
x -> y [1].

# References
1. Paluš, M., Krakovská, A., Jakubík, J., & Chvosteková, M. (2018). Causality,
dynamical systems and the arrow of time. Chaos: An Interdisciplinary Journal of
Nonlinear Science, 28(7), 075307. http://doi.org/10.1063/1.5019944
"""
function eom_ar1(x, p, n)
    a₁, b₁, c₁, σ = (p...,)
    x, y = (x...,)
    ξ₁ = rand(Normal(0, σ))
    ξ₂ = rand(Normal(0, σ))

    dx = a₁*x + ξ₁
    dy = b₁*y + c₁*x + ξ₂
    return SVector{2}(dx, dy)
end

function ar1(uᵢ, a₁, b₁, c₁, σ)
    p = [a₁, b₁, c₁, σ]
    return DiscreteDynamicalSystem(eom_ar1, uᵢ, p)
end

"""
    ar1(uᵢ, a₁ = 0.90693, b₁ = 0.40693, c₁ = 0.5, σ = 0.40662) -> DiscreteDynamicalSystem

Initialise a discrete vector autoregressive system where x → y
"""
ar1(;uᵢ = rand(2), a₁ = 0.90693, b₁ = 0.40693, c₁ = 0.5, σ = 0.40662) =
    ar1(uᵢ, a₁, b₁, c₁, σ)
