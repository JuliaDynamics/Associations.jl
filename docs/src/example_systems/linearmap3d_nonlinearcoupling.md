# Synthetic coupled dynamical systems

## Three linear maps with nonlinear coupling

For this example, we'll consider a 3d linear system (``x_1``, ``x_2`` and ``x_3``) with nonlinear coupling,  where the dynamical influence goes in the directions ``x_1 \rightarrow x_2``, ``x_1 \rightarrow x_3`` and ``x_2 \rightarrow x_3``.

The system is given by the following difference equations:

```math
\small
\begin{aligned}
x_1(t+1) &= a_1 x_1 (1-x_1(t))^2  e^{-x_2(t)^2} + 0.4 \xi_{1}(t) \\
x_2(t+1) &= a_1 x_2 (1-x_2(t))^2  e^{-x_2(t)^2} + 0.4 \xi_{2}(t) + b x_1 x_2 \\
x_3(t+1) &= a_3 x_3 (1-x_3(t))^2  e^{-x_3(t)^2} + 0.4 \xi_{3}(t) + c x_{2}(t) + d x_{1}(t)^2.
\end{aligned}
\normalsize
```

Here, ``\xi_{1,2,3}(t)`` are independent normally distributed noise processes, representing dynamical noise in the system, with zero mean and standard deviations ``\sigma_1``, ``\sigma_2``, ``\sigma_3``,
respectively.

## Where has the system been used?
This system was used in [Gourévitch et al. (2006)](https://link.springer.com/article/10.1007/s00422-006-0098-0),
 where they review the linear Granger causality test, partial directed coherence and the nonlinear Granger causality test.

## Represent as a DiscreteDynamicalSystem

```@setup linear3d_nonlinearcoupling
using CausalityTools
using DynamicalSystems
using Plots
using StaticArrays
using Statistics
using Distributions
```

We first define the equations of motion.

```@example linear3d_nonlinearcoupling
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
nothing #hide
```

To make things easier to use, we create function that generates a
DiscreteDynamicalSystem instance for any set of parameters `a₁`, `a₂`, `a₃`, `b`, `c` and `d`, initial condition `u₀`, and dynamical noise levels `σ₁`, `σ₂` and `σ₃`.

```@example linear3d_nonlinearcoupling
function linear3d_nonlinearcoupling(;uᵢ = rand(3),
                σ₁ = 1.0, σ₂ = 1.0, σ₃ = 1.0,
                a₁ = 3.4, a₂ = 3.4, a₃ = 3.4,
                b = 0.5, c = 0.3, d = 0.5)
    p = [a₁, a₂, a₃, b, c, d, σ₁, σ₂, σ₃]
    return DiscreteDynamicalSystem(eom_linear3d_nonlinearcoupling, uᵢ, p)
end
nothing #hide
```

An example realization of the system is:

```@example linear3d_nonlinearcoupling
s = linear3d_nonlinearcoupling()
orbit = trajectory(s, 100)
x1, x2, x3 = orbit[:, 1], orbit[:, 2], orbit[:, 3]
plot(x1, label = "x", lc = :black)
plot!(x2, label = "y", lc = :red)
plot!(x3, label = "z", lc = :blue)
xlabel!("Time step"); ylabel!("Value")
savefig("linear3d_nonlinearcoupling.svg") #hide
nothing #hide
```

![](linear3d_nonlinearcoupling.svg)


## Predefined system
This system is predefined in `CausalityTools.Systems`, and can be initialized using the `linearmap3d_nonlinearcoupling` function.


## References
Gourévitch, B., Le Bouquin-Jeannès, R., & Faucon, G. (2006). Linear and
nonlinear causality between signals: methods, examples and neurophysiological applications. Biological Cybernetics, 95(4), 349–369.
[https://link.springer.com/article/10.1007/s00422-006-0098-0](https://link.springer.com/article/10.1007/s00422-006-0098-0)
