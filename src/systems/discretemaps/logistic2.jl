"""
    eom_logistic2(dx, x, p, n) -> function

Equations of motions for a system consisting of two coupled logistic maps where
X unidirectionally influences Y.

The parameter `c` controls how strong the dynamical forcing is. Parameters `r₁`
and `r₂` are set to the chaotic regime by default

# References
1. Diego, D., Agasøster Haaga, K., & Hannisdal, B. (2018, November 1).
Transfer entropy computation using the Perron-Frobenius operator.
Eprint ArXiv:1811.01677. Retrieved from
https://ui.adsabs.harvard.edu/#abs/2018arXiv181101677D
"""
function eom_logistic2(dx, x, p, n)
    c, r₁, r₂, σ = (p...,)
    ξ = rand() # random number from flat distribution on [0, 1]
    x, y = x[1], x[2]
    f_xy = (y +  (c*(x + σ*ξ)/2) ) / (1 + (c/2)*(1+σ))

    dx[1] = r₁ * x * (1 - x)
    dx[2] = r₂ * (f_xy) * (1 - f_xy)
    return
end

function logistic2(u₀, c, r₁, r₂, σ)
    p = [c, r₁, r₂, σ]
    DiscreteDynamicalSystem(eom_logistic2, u₀, p)
end

"""
    logistic2(;u₀ = rand(2), c = 0.1, σ = 0.05,
        r₁ = 3.78, r₂ = 3.66) -> DiscreteDynamicalSystem

Initialise a system consisting of two coupled logistic maps where X
unidirectionally influences Y. By default, the parameters `r₁` and `r₂` are set
to values yielding chaotic behaviour.

The parameter `c` controls how strong the dynamical forcing is. If `σ > 0`,
dynamical noise masking the influence of  `x` on `y` equivalent to
``\\sigma \\cdot \\xi`` is added at each iteration. Here,``\\xi`` is a draw from a
flat distribution on ``[0, 1]``. Thus, setting `σ = 0.05` is equivalent to
add dynamical noise corresponding to a maximum of ``5 \\%`` of the possible
range of values of the logistic map.



The equations of motion are

```math
\\begin{aligned}
dx &= r_1 x(1 - x) \\
dy &= r_2 f(x,y)(1 - f(x,y)),
\\end{aligned}
```
with
```math
\\begin{aligned}
f(x,y) = \\dfrac{y + \\frac{c(x \\xi )}{2}}{1 + \\frac{c}{2}(1+ \\sigma )}
\\end{aligned}
```

where

# References

1. Diego, D., Agasøster Haaga, K., & Hannisdal, B. (2018, November 1).
Transfer entropy computation using the Perron-Frobenius operator.
Eprint ArXiv:1811.01677. Retrieved from
https://ui.adsabs.harvard.edu/#abs/2018arXiv181101677D

"""
logistic2(;u₀ = rand(2), c = 0.1, r₁ = 3.78, r₂ = 3.66, σ = 0.05) =
    logistic2(u₀, c, r₁, r₂, σ)
