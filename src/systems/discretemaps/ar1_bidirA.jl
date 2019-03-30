"""
    eom_ar1_bidirA(x, p, n) -> SVector{2}

Equations of motion a system consisting of two mutually coupled first order
autoregressive processes. For `c = 0`, only Y -> X. For `c > 0`, there is
bidirectional coupling. Each subcomponent is dynamically affected by
an independent Gaussian random processes (with standard deviations `σx` and
`σy`) at each time step.

"""
function eom_ar1_bidirA(x, p, n)
    c_xy, c_yx, σx, σy = (p...,)
    x, y = (x...,)

    ϵx = rand(Normal(0, σx))
    ϵy = rand(Normal(0, σy))

    dx = 0.5*x + c_yx*y + ϵx
    dy = c_xy*x + 0.7*y + ϵy
    return SVector{2}(dx, dy)
end

function ar1_bidirA(uᵢ, c_xy, c_yx, σx, σy)
    p = [c_xy, c_yx, σx, σy]
    logistic_system = DiscreteDynamicalSystem(eom_ar1_bidirA, uᵢ, p)

    return logistic_system
end

"""
    ar1_bidirA(;uᵢ = rand(2), c_xy = 0, c_yx = 0.2, σx = sqrt(0.1), σy = sqrt(0.1))

A system consisting of two mutually coupled first order autoregressive processes. 
For `c = 0`, only Y -> X. For `c > 0`, there is bidirectional coupling. Each subcomponent 
is dynamically affected by an independent Gaussian random processes (with standard 
deviations `σx` and `σy`) at each time step.
"""
ar1_bidirA(;uᵢ = rand(2), c_xy = 0, c_yx = 0.2, σx = sqrt(0.1), σy = sqrt(0.1)) =
    ar1_bidirA(uᵢ, c_xy, c_yx, σx, σy)
