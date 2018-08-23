"""
    _logistic2(dx, x, p, n) -> function

Equations of motion for a unidirectionally coupled logistic map system where
`x₁` drives `x₂`. This example is from [1].

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

"""
     logistic2(u₀, c)

Initialise a unidirectionally coupled logistic map system where `x₁` drives
`x₂`. This example is from [1].

# References
1. D Diego, KA Haaga, B Hannisdal, in prep. Transfer Entropy computation by
Perron-Frobenius operator approximation.
"""
function logistic2(u₀, c, r₁, r₂)
    p = [c, r₁, r₂]
    DiscreteDynamicalSystem(eom_logistic2, u₀, p)
end

logistic2(;u₀ = rand(2), c = 2.0, r₁ = 3.78, r₂ = 3.66) =
    logistic2(u₀, c, r₁, r₂)
