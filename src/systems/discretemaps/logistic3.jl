
doc"""
    eom_logistic3(u, p, t)

Equations of motion for a triple coupled logistic map system with a common
driver, where `Z` -> `X` and `Z -> Y`. This is equation 36 in [1].


# References
1. Runge, Jakob. Causal network reconstruction from time series: From
theoretical assumptions to practical estimation, Chaos 28, 075310 (2018);
doi: 10.1063/1.5025050
"""
function eom_logistic3(u, p, t)
    r, σx, σy, σz = (p...)
    x, y, z = (u...)

    # Independent dynamical noise for each variable.
    ηx = rand()
    ηy = rand()
    ηz = rand()

    dx = (x*(r - r*x - z + σx*ηx)) % 1
    dy = (y*(r - r*y - z + σy*ηy)) % 1
    dz = (z*(r - r*z + σz*ηz)) % 1
    return SVector{3}(dx, dy, dz)
end


function logistic3(u₀, r, σx, σy, σz)
    p = [r, σx, σy, σz]
    DiscreteDynamicalSystem(eom_logistic3, u₀, p)
end

"""
    logistic3(;u₀ = rand(3), r = 4, σx = 0.05, σy = 0.05, σz = 0.05)

Initialise a triple coupled logistic map system with a common driver, where
`Z` -> `X` and `Z -> Y` [1].

Dynamical noise may be added to each of the dynamical variables by tuning the
parameters `σz`, `σx` and `σz`. The parameter `r` is set to 4 by default,
which is in the chaotic regime.

# References
1. Runge, Jakob. Causal network reconstruction from time series: From
theoretical assumptions to practical estimation, Chaos 28, 075310 (2018);
doi: 10.1063/1.5025050
"""
logistic3(;u₀ = rand(3), r = 4, σx = 0.05, σy = 0.05, σz = 0.05) =
    logistic3(u₀, r, σx, σy, σz)
