function eom_nonlinear_sindriver_2d(dx, x, p, n)
    a, b, c, t, Δt = (p...,)
    x, y = x[1], x[2]
    𝒩 = Normal(0, 1)
    
    dx[1] = sin(t)
    dx[2] = a*x * (1 - b*x) + c*rand(𝒩)
    p[end-1] += Δt # update t

    return
end

"""
    nonlinear_sindriver_2d(;u₀ = 0, a = 1.0, b = 1.0, c = 2.0, Δt = 1)

A discrete, nonlinear 2D system from McCracken & Weigel (2016)[^McCrackenWeigel2016],
governed by the following equations of motion:

```math
\\begin{aligned}
x_{t+1} = \\sin(t) \\\\
y(t+1) &= a x_t (1 - B x_t) + C \\eta_t
\\end{aligned}
```

where ``\\eta_t \\sim \\mathcal{N}(\\mu = 0, \\sigma = 1)``, ``A, B, C \\in [0, 1]``.

[^McCrackenWeigel2016]: McCracken, J. M., & Weigel, R. S. (2016). Nonparametric causal inference for bivariate time series. Physical Review E, 93(2), 022207.
"""
function nonlinear_sindriver_2d(;u₀ = 0, a = 1.0, b = 1.0, c = 2.0, Δt = (1/30*π))
    @assert a >= 0; @assert b >= 0; @assert c >= 0

    DiscreteDynamicalSystem(eom_nonlinear_sindriver2d, [0, u₀], [a, b, c, 0, Δt])
end