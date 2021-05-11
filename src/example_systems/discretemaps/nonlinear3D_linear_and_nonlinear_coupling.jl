using LabelledArrays

export nonlinear3d

"""
    eom_nonlinear3d(u₀, a₁, a₂, a₃,  b₁, b₂, b₃, 
        c₁₂, c₂₃, c₁₃, σ₁, σ₂, σ₃) → DiscreteDynamicalSystem

Equations of motion for a 3d nonlinear system with nonlinear couplings 
``x_1 \\to x_2``, ``x_2 \\to x_3`` and ``x_1 \\to x_3``. Modified from [1]. 

## Equations of motion 

The equations of motion are

```math
\\begin{aligned}
x_1(t+1) &= a_1 x_1 (1-x_1(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{1}(t) \\\\
x_2(t+1) &= a_1 x_2 (1-x_2(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{2}(t) + b x_1 x_2 \\\\
x_3(t+1) &= a_3 x_3 (1-x_3(t))^2  e^{-x_3(t)^2} + 0.4 \\xi_{3}(t) + c x_{2}(t) + d x_{1}(t)^2.
\\end{aligned}
```

## References 

1. Gourévitch, B., Le Bouquin-Jeannès, R., & Faucon, G. (2006). Linear and nonlinear 
    causality between signals: methods, examples and neurophysiological 
    applications. Biological Cybernetics, 95(4), 349–369.
"""
function eom_nonlinear3d(x, p, n)
    x₁, x₂, x₃ = (x...,)
    a₁, a₂, a₃, b₁, b₂, b₃, c₁₂, c₂₃, c₁₃, σ₁, σ₂, σ₃ = (p...,)
    ξ₁ = rand(Normal(0, σ₁))
    ξ₂ = rand(Normal(0, σ₂))
    ξ₃ = rand(Normal(0, σ₃))
    
    dx₁ = a₁*x₁*(1-x₁)^2 * exp(-x₁^2) + b₁*ξ₁
    dx₂ = a₂*x₂*(1-x₂)^2 * exp(-x₂^2) + b₂*ξ₂ + c₁₂*x₁*x₂ 
    dx₃ = a₃*x₃*(1-x₃)^2 * exp(-x₃^2) + b₃*ξ₃ + c₂₃*x₂ + c₁₃*x₁^2

    return SVector{3}(dx₁, dx₂, dx₃)
end

function nonlinear3d(u₀, a₁, a₂, a₃,  b₁, b₂, b₃, c₁₂, c₂₃, c₁₃, σ₁, σ₂, σ₃)
    p = @LArray [a₁, a₂, a₃,  b₁, b₂, b₃, c₁₂, c₂₃, c₁₃, σ₁, σ₂, σ₃] (:a₁, :a₂, :a₃,  :b₁, :b₂, :b₃, :c₁₂, :c₂₃, :c₁₃, :σ₁, :σ₂, :σ₃)
    s = DiscreteDynamicalSystem(eom_nonlinear3d, u₀, p)
    return s
end

"""
    nonlinear3d(;u₀ = rand(3), 
        σ₁ = 1.0, σ₂ = 1.0, σ₃ = 1.0, 
        a₁ = 3.4, a₂ = 3.4, a₃ = 3.4, 
        b₁ = 0.4, b₂ = 0.4, b₃ = 0.4, 
        c₁₂ = 0.5, c₂₃ = 0.3, c₁₃ = 0.5) → DiscreteDynamicalSystem

A 3d nonlinear system with nonlinear couplings ``x_1 \\to x_2``, 
``x_2 \\to x_3`` and ``x_1 \\to x_3``. Modified from [1]. 

## Equations of motion 

The equations of motion are

```math
\\begin{aligned}
x_1(t+1) &= a_1 x_1 (1-x_1(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{1}(t) \\\\
x_2(t+1) &= a_1 x_2 (1-x_2(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{2}(t) + b x_1 x_2 \\\\
x_3(t+1) &= a_3 x_3 (1-x_3(t))^2  e^{-x_3(t)^2} + 0.4 \\xi_{3}(t) + c x_{2}(t) + d x_{1}(t)^2.
\\end{aligned}
```

## References 

1. Gourévitch, B., Le Bouquin-Jeannès, R., & Faucon, G. (2006). Linear and nonlinear 
    causality between signals: methods, examples and neurophysiological 
    applications. Biological Cybernetics, 95(4), 349–369.
"""
nonlinear3d(;u₀ = rand(3), 
        σ₁ = 1.0, σ₂ = 1.0, σ₃ = 1.0, 
        a₁ = 3.4, a₂ = 3.4, a₃ = 3.4, 
        b₁ = 0.4, b₂ = 0.4, b₃ = 0.4, 
        c₁₂ = 0.5, c₂₃ = 0.3, c₁₃ = 0.5) = 
    nonlinear3d(u₀, a₁, a₂, a₃,  b₁, b₂, b₃, c₁₂, c₂₃, c₁₃, σ₁, σ₂, σ₃)
