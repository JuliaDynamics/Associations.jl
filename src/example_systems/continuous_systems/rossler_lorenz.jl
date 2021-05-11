using LabelledArrays

@inline @inbounds function eom_rossler_lorenz(u, p, t)
    c_xy, a₁, a₂, a₃, b₁, b₂, b₃ = (p...,)
    x1, x2, x3, y1, y2, y3 = u[1], u[2], u[3], u[4], u[5], u[6] 
    
    dx1 = -a₁*(x2 + x3)
    dx2 = a₁*(x1 + a₂*x2)
    dx3 = a₁*(a₂ + x3*(x1 - a₃))
    dy1 = b₁*(-y1 + y2)
    dy2 = b₂*y1 - y2 - y1*y3 + c_xy*(x2^2)
    dy3 = y1*y2 - b₃*y3
    
    return SVector{6}(dx1, dx2, dx3, dy1, dy2, dy3)
end

"""
    rossler_lorenz(;u₀ = rand(6), a₁ = -6, a₂ = 6, a₃ = 2.0, 
        b₁ = 10, b₂ = 28, b₃ = 8/3, c_xy = 1) → ContinuousDynamicalSystem

Initialise a Rössler-Lorenz system consisting of two independent 3D subsystems:
one Rössler system and one Lorenz system. They are coupled such that the
second component (`x₂`) of the Rössler system unidirectionally forces the
second component (`y₂`) of the Lorenz system. 

The parameter `c_xy` controls the coupling strength. The implementation here also 
allows for tuning the parameters of each subsystem by introducing the constants 
`a₁`, `a₂`, `a₃`, `b₁`, `b₂`, `b₃`. Default values for these parameters are 
as in [1].

## Equations of motion 

The dynamics is generated by the following vector field

```math
\\begin{aligned}
\\dot x_1 &= a_1(x_2 + x_3) \\\\
\\dot x_2 &= a_2(x_1 + 0.2x_2) \\\\
\\dot x_3 &= a_2(0.2 + x_3(x_1 - a_3)) \\\\
\\dot y_1 &= b_1(y_2 - y_1) \\\\
\\dot y_2 &= y_1(b_2 - y_3) - y_2 +c_{xy}(x_2)^2 \\\\
\\dot y_3 &= y_1 y_2 - b_3y_3
\\end{aligned}
```

with the coupling constant ``c_{xy} \\geq 0``.

## References

1. Krakovská, Anna, et al. "Comparison of six methods for the detection of causality in a 
    bivariate time series." Physical Review E 97.4 (2018):042207. 
    [https://journals.aps.org/pre/abstract/10.1103/PhysRevE.97.042207](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.97.042207)
"""
function rossler_lorenz(;u₀ = rand(6), a₁ = -6, a₂ = 6, a₃ = 2.0, 
    b₁ = 10, b₂ = 28, b₃ = 8/3, c_xy = 1)

    p = @LArray [c_xy, a₁, a₂, a₃, b₁, b₂, b₃] (:c_xy, :a₁, :a₂, :a₃, :b₁, :b₂, :b₃)
    ContinuousDynamicalSystem(eom_rossler_lorenz, u₀, p)
end
