doc"""
    _logistic2(dx, x, p, n) -> function

Equations of motions for a system consisting of two coupled logistic maps where
X unidirectionally influences Y.

The parameter `c` controls how strong the dynamical forcing is. Parameters `r₁`
and `r₂` are set to the chaotic regime by default

# References
1. D Diego, KA Haaga, B Hannisdal, in prep. Transfer Entropy computation by
Perron-Frobenius operator approximation.
"""
function eom_logistic2(dx, x, p, n)
    c, r₁, r₂ = (p...)
    x, y = x[1], x[2]
    f_xy = (y +  (c * x/2)) / (1 + (c/2))

    dx[1] = r₁ * x * (1 - x)
    dx[2] = r₂ * (f_xy) * (1 - f_xy)
    return
end

function logistic2(u₀, c, r₁, r₂)
    p = [c, r₁, r₂]
    DiscreteDynamicalSystem(eom_logistic2, u₀, p)
end

doc"""
    logistic2(;u₀ = rand(2), c = 2.0,
        r₁ = 3.78, r₂ = 3.66) -> DiscreteDynamicalSystem

Initialise a system consisting of two coupled logistic maps where X
unidirectionally influences Y.

The parameter `c` controls how strong the dynamical forcing is. Parameters `r₁`
and `r₂` are set to the chaotic regime by default.

The equations of motion are

```math
\begin{aligned}
dx &= r_1x(1 - x)
dy &= r_2f(x,y)(1 - f(x,y)),
\end{aligned}
```
with
```math
\begin{aligned}
f(x,y) = \dfrac{y + \frac{cx}{2}}{1 + \frac{c}{2}}
\end{aligned}
```

# References
D Diego, KA Haaga, B Hannisdal, in prep. Transfer Entropy computation by Perron-Frobenius operator approximation.
"""
logistic2(;u₀ = rand(2), c = 2.0, r₁ = 3.78, r₂ = 3.66) =
    logistic2(u₀, c, r₁, r₂)
