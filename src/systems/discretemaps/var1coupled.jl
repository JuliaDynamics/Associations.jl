"""
    eom_var1coupled(x, p, n) -> SVector{3}

Equations of motion a system consisting of two mutually coupled first order
autoregressive processes. For `c = 0`, only Y -> X. For `c > 0`, there is
bidirectional coupling. Each subcomponent is dynamically affected by
an independent Gaussian random processes (with standard deviations `σx` and
`σy`) at each time step.

"""
function eom_var1coupled(x, p, n)
    c, σx, σy = (p...,)
    x, y = (x...,)

    ϵx = rand(Normal(0, σx))
    ϵy = rand(Normal(0, σy))

    dx = 0.5*x + 0.2*y + ϵx
    dy = c*x + 0.7*y + ϵy
    return SVector{2}(dx, dy)
end

function var1coupled(uᵢ, c, σx, σy)
    p = [c, σx, σy]
    logistic_system = DiscreteDynamicalSystem(eom_var1coupled, uᵢ, p)

    return logistic_system
end

"""
    var1coupled(;uᵢ = rand(2), c = 0, σx = sqrt(0.1), σy = sqrt(0.1))

Initialise a system consisting of two mutually coupled first order
autoregressive processes. For `c = 0`, only Y -> X. For `c > 0`, there is
bidirectional coupling. Each subcomponent is dynamically affected by
an independent Gaussian random processes (with standard deviations `σx` and
`σy`) at each time step.
"""
var1coupled(;uᵢ = rand(2), c = 0, σx = sqrt(0.1), σy = sqrt(0.1)) =
var1coupled(uᵢ, c, σx, σy)
