"""
    eom_linear3d_nonlinearcoupling(x, p, n) -> Function


Equations of motion for a 3d linear system with nonlinear coupling [1].
The difference equations are

```math
\\begin{aligned}
x_1(t+1) &= a_1 x_1 (1-x_1(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{1}(t) \\
x_2(t+1) &= a_1 x_2 (1-x_2(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{2}(t) + b x_1 x_2 \\
x_3(t+1) &= a_3 x_3 (1-x_3(t))^2  e^{-x_3(t)^2} + 0.4 \\xi_{3}(t) + c x_{2}(t) \\
         &+ d x_{1}(t)^2.
\\end{aligned}
```

Here, ``\\xi_{1,2,3}(t)`` are independent normally distributed noise processes,
representing dynamical noise in the system, with zero mean and standard
deviations ``\\sigma_1``, ``\\sigma_2``, ``\\sigma_3``, respectively.


# References
1. Gourévitch, B., Le Bouquin-Jeannès, R., & Faucon, G. (2006). Linear and
nonlinear causality between signals: methods, examples and neurophysiological
applications. Biological Cybernetics, 95(4), 349–369.
"""
function eom_linear3d_nonlinearcoupling(x, p, n)
    x₁, x₂, x₃ = (x...,)
    a₁, a₂, a₃, b, c, d, σ₁, σ₂, σ₃ = (p...,)
    ξ₁ = rand(Normal(0, σ₁))
    ξ₂ = rand(Normal(0, σ₂))
    ξ₃ = rand(Normal(0, σ₃))

    dx₁ = a₁*x₁*(1-x₁)^2 * exp(-x₁^2) + 0.4*ξ₁
    dx₂ = a₂*x₂*(1-x₂)^2 * exp(-x₂^2) + 0.4*ξ₂ + b*x₁*x₂
    dx₃ = a₃*x₃*(1-x₃)^2 * exp(-x₃^2) + 0.4*ξ₃ + c*x₂ + d*x₁^2

    return SVector{3}(dx₁, dx₂, dx₃)
end


function linear3d_nonlinearcoupling(uᵢ, a₁, a₂, a₃, b, c, d, σ₁, σ₂, σ₃)
    p = [a₁, a₂, a₃, b, c, d, σ₁, σ₂, σ₃]
    s = DiscreteDynamicalSystem(eom_linear3d_nonlinearcoupling, uᵢ, p)
    return s
end

"""
    linear3d_nonlinearcoupling(;uᵢ = rand(3), σ₁ = 1.0, σ₂ = 1.0, σ₃ = 1.0,
        a₁ = 3.4, a₂ = 3.4, a₃ = 3.4, b = 0.5, c = 0.3, d = 0.5) -> DiscreteDynamicalSystem

A 3d linear system with nonlinear coupling [1]. The difference equations are

```math
\begin{aligned}
x_1(t+1) &= a_1 x_1 (1-x_1(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{1}(t) \\
x_2(t+1) &= a_1 x_2 (1-x_2(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{2}(t) + b x_1 x_2 \\
x_3(t+1) &= a_3 x_3 (1-x_3(t))^2  e^{-x_3(t)^2} + 0.4 \\xi_{3}(t) + c x_{2}(t) \\
         &+ d x_{1}(t)^2.
\end{aligned}
```

Here, ``\\xi_{1,2,3}(t)`` are independent normally distributed noise processes,
representing dynamical noise in the system, with zero mean and standard
deviations ``\\sigma_1``, ``\\sigma_2``, ``\\sigma_3``, respectively.


# References
1. Gourévitch, B., Le Bouquin-Jeannès, R., & Faucon, G. (2006). Linear and
nonlinear causality between signals: methods, examples and neurophysiological
applications. Biological Cybernetics, 95(4), 349–369.
"""
linear3d_nonlinearcoupling(;uᵢ = rand(3), σ₁ = 1.0, σ₂ = 1.0, σ₃ = 1.0,
        a₁ = 3.4, a₂ = 3.4, a₃ = 3.4, b = 0.5, c = 0.3, d = 0.5) =
    linear3d_nonlinearcoupling(uᵢ, a₁, a₂, a₃, b, c, d, σ₁, σ₂, σ₃)
