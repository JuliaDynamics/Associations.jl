using LabelledArrays

export ar1_unidir

"""
    eom_ar1(x, p, n) → Function

Equations of motion for a bivariate, order one autoregressive model, where
x → y [1].

## Equations of motion

```math
\\begin{aligned}
x(t+1) &= a_1 x(t) + \\xi_{1} \\\\
y(t+1) &= b_1 y(t) - c_{xy} x + \\xi_{2},
\\end{aligned}
```

where ``\\xi_{1}`` and ``\\xi_{2}`` are drawn from normal distributions 
with zero mean and standard deviation `σ` at each iteration.

## References

1. Paluš, M., Krakovská, A., Jakubík, J., & Chvosteková, M. (2018). Causality,
    dynamical systems and the arrow of time. Chaos: An Interdisciplinary Journal of
    Nonlinear Science, 28(7), 075307. http://doi.org/10.1063/1.5019944
"""
function eom_ar1_unidir(x, p, n)
    a₁, b₁, c_xy, σ = (p...,)
    x, y = (x...,)
    ξ₁ = rand(Normal(0, σ))
    ξ₂ = rand(Normal(0, σ))
    
    dx = a₁*x + ξ₁
    dy = b₁*y + c_xy*x + ξ₂
    return SVector{2}(dx, dy)
end

function ar1_unidir(u₀, a₁, b₁, c_xy, σ)
    p = @LArray [a₁, b₁, c_xy, σ] (:a₁, :b₁, :c_xy, :σ)
    DiscreteDynamicalSystem(eom_ar1_unidir, u₀, p)
end

"""
    ar1_unidir(u₀, a₁ = 0.90693, b₁ = 0.40693, c_xy = 0.5, 
        σ = 0.40662) → DiscreteDynamicalSystem

A bivariate, order one autoregressive model, where ``x \\to y`` [1].

## Equations of motion

```math
\\begin{aligned}
x(t+1) &= a_1 x(t) + \\xi_{1} \\\\
y(t+1) &= b_1 y(t) - c_{xy} x + \\xi_{2},
\\end{aligned}
```

where ``\\xi_{1}`` and ``\\xi_{2}`` are drawn from normal distributions 
with zero mean and standard deviation `σ` at each iteration.

## References

1. Paluš, M., Krakovská, A., Jakubík, J., & Chvosteková, M. (2018). Causality,
    dynamical systems and the arrow of time. Chaos: An Interdisciplinary Journal of
    Nonlinear Science, 28(7), 075307. http://doi.org/10.1063/1.5019944
"""
ar1_unidir(;u₀ = rand(2), a₁ = 0.90693, b₁ = 0.40693, c_xy = 0.5, σ = 0.40662) = 
    ar1_unidir(u₀, a₁, b₁, c_xy, σ)