"""
    eom_henon4(x, p, n) -> Function

Equations of motion for a 4D Henon system consisting of two identical
Henon maps, `X` and `Y` with unidirectional forcing from `X` to`Y` [1].
The coupling constant `c` controls the strength of the forcing.

The implementation here has adjustable parameters; the default values are the
ones used in [1].

# References
1. Krakovská, A., Jakubík, J., Chvosteková, M., Coufal, D.,
Jajcay, N., & Paluš, M. (2018). Comparison of six methods for the detection of
causality in a bivariate time series. Physical Review E, 97(4), 042207.
"""
function eom_henon4(x, p, n)
    a₁, a₂, b₁, b₂, c, k = (p...,)
    x₁, x₂, y₁, y₂ = (x...,)

    dx₁ = a₁ - x₁^2 + k*x₂
    dx₂ = b₁*x₁
    dy₁ = a₂ - (c * x₁ * y₁ + (1 - c)*y₁^2) + k*y₂
    dy₂ = b₂*y₁
    return SVector{4}(dx₁, dx₂, dy₁, dy₂)
end

"""
    henon4(u₀, c) -> DiscreteDynamicalSystem

Initialize an instance of a 4D Henon map system consisting of two identical
Henon systems, ``x`` and ``y`` [1]. The subsystems are unidirectionally
coupled . The coupling constant `c` controls the strength of the forcing.
Synchronization occurs when the values of the coupling constant
``c > 0.7``.

The difference equations are:

```math
\\begin{aligned}
x_1(t+1) &= a_1 - x_1^2(t) + k x_2(t) \\
x_2(t+1) &= x_1(t) \\
y_1(t+1) &= a_1 - [c x_1(t) y_1(t) + (1-c) y_1^2(t)] + k y_2(t) \\
y_2(t+1) &= y_1(t)
\\end{aligned}
```

This system was investigated by Krakovská to study the performance of different
causality detection algorithms.

# References
Krakovská, A., Jakubík, J., Chvosteková, M., Coufal, D., Jajcay, N., & Paluš, M. (2018). Comparison of six methods for the detection of causality in a bivariate time series. Physical Review E, 97(4), 042207.
"""
function henon4(u₀, a₁, a₂, b₁, b₂, c, k)
    p = [a₁, a₂, b₁, b₂, c, k]
    logistic_system = DiscreteDynamicalSystem(eom_henon4, u₀, p)
    return logistic_system
end

henon4(;u₀ = rand(4),
        a₁ = 1.4, a₂ = 1.4, b₁ = 1.0, b₂ = 1.0, c = 2.0, k = 0.3) =
    henon4(u₀, a₁, a₂, b₁, b₂, c, k)
