export lorenzdiffusive

import StaticArrays: SVector
import DynamicalSystemsBase: ContinuousDynamicalSystem

export lorenzdiffusive

@inline @inbounds function eom_lorenzdiffusive(u, p, t)
    
    C₁₂, C₂₁, R, ϵ₁, ϵ₂ = p[1], p[2], p[3], p[4], p[5]
    x₁, x₂, x₃ = u[1], u[2], u[3]
    y₁, y₂, y₃ = u[4], u[5], u[6]

    dx1 = 10*(x₂ - x₁) + C₂₁*(y₁-x₁)
    dx2 = x₁*((R+ϵ₁) - x₃) - x₂
    dx3 = x₁*x₂ - 8/3*x₃
    
    dy1 = 10*(y₂ - y₁) + C₁₂*(x₁-y₁)
    dy2 = y₁*((R+ϵ₂) - y₃) - y₂
    dy3 = y₁*y₂ - 8/3*y₃
    return SVector{6}(dx1, dx2, dx3, dy1, dy2, dy3)
end

"""
    lorenzdiffusive(; ui = rand(6), C₁₂::Real = 5, C₂₁::Real = 0, 
        R::Real = 28.0, ϵ₁::Real = -0.02, ϵ₂::Real = 0.03)

A dynamical system consisting of two diffusively coupled 3D Lorenz systems[^Martini2011].

The coupling magnitude from subsystem 1 to subsystem 2 is controlled by `C₁₂`, and the 
coupling from subsystem 2 to subsystem 1 is controlled by `C₂₁`. The parameters `ϵ₁` and `ϵ₂` 
add small deviations to the control parameter `R`.

## Equations of motion

```math 
\\begin{aligned}
\\dot{x_1} &= 10(x_2 - x_1) + C_{21}*(y_1-x_1) \\\\
\\dot{x_2} &= x_1((R+ϵ₁) - x_3) - x_2 \\\\
\\dot{x_3} &= x_1x_2 - 8/3x_3 \\\\
\\dot{y_1} &= 10(y_2 - y_1) + C_{12}(x_1-y_1) \\\\
\\dot{y_2} &= y_1((R+\\epsilon_2) - y_3) - y_2 \\\\
\\dot{y_3} &= y_1y_2 - 8/3y_3
\\end{aligned}
```

[^Martini2011]: Martini, M., Kranz, T. A., Wagner, T., & Lehnertz, K. (2011). Inferring directional interactions from transient signals with symbolic transfer entropy. Physical review E, 83(1), 011919.
"""
function lorenzdiffusive(; ui = rand(6), C₁₂::Real = 5, C₂₁::Real = 0, 
        R::Real = 28.0, ϵ₁::Real = -0.02, ϵ₂::Real = 0.03)

    p = [C₁₂, C₂₁, R, ϵ₁, ϵ₂]
    ContinuousDynamicalSystem(eom_lorenzdiffusive, ui, p)
end

