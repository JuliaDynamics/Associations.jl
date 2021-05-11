export hindmarsh_rose

"""
Equations of motion for the Hindmarsh-Rose system.
"""
function eom_hindmarsh_rose(u, p, t)
    a, b, c, d, r, s, xᵣ, I = (p...,)
    x, y, z = (u...,)

	ϕ = -a*x^3 + b*x^2
	ψ = c - d*x^2
    dx = y + ϕ - z + I
	dy = ψ - y
	dz = r*(s*(x - xᵣ) - z)
    return SVector{3}(dx, dy, dz)
end

"""
	hindmarsh_rose(u₀, p)

Initialise a Hindmarsh-Rose system, which is a model of neuronal
spiking.

```math
\\begin{aligned}
\\dfrac{dx}{dt} &= y + \\phi(x) - z + I
\\dfrac{dy}{dt} &= \\psi(x) - y
\\dfrac{dz}{dt} &= r[s(x - x_R) - z],
\\end{aligned}
```
where

```math
\\begin{aligned}
\\phi(x) &= -ax^3+bx^2
\\psi(x) &= c - dx^2
\\end{aligned}
```

If parameters other than the defaults are to be used, they must be
provided as a vector `[a, b, c, d, r, s, xᵣ, I]`.
"""
function hindmarsh_rose(u₀, p)
    ContinuousDynamicalSystem(eom_hindmarsh_rose, u₀, p)
end
hindmarsh_rose(;u₀ = rand(3), a = 1, b = 3, c = 1, d = 5, r = 1e-3, s = 4, xᵣ = - 8/5, I = -8) =
    hindmarsh_rose(u₀, [a, b, c, d, r, s, xᵣ, I])
