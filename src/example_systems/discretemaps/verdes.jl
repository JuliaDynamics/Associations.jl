using LabelledArrays

export verdes

"""
    eom_verdes(u, p, t) → SVector{3}

Equations of motion for a 3D system where the response X is a highly
nonlinear combination of Y and Z [1]. The forcings Y and Z involve
sines and cosines, respectively, and have different periods.

The implementation here allows for tuning the period of the
forcing signals.

# References
Verdes, P. F. "Assessing causality from multivariate time series." Physical Review E 72.2 (2005): 026222.
"""
function eom_verdes(u, p, t)
    x, y, z = (u...,)
    ωy, ωz, σx, σy, σz = (p...,)

    ηx = σx == 0 ? 0 : rand(Normal(0, σx))
    ηy = σy == 0 ? 0 : rand(Normal(0, σy))
    ηz = σz == 0 ? 0 : rand(Normal(0, σz))

    dx = y*(18y - 27y^2 + 10)/2 + z*(1-z) + ηx
    dy = (1 - cos((2*pi/ωy) * t))/2 + ηy
    dz = (1 - sin((2*pi/ωz) * t))/2 + ηz
    return SVector{3}(dx, dy, dz)
end


function verdes(u₀, ωy, ωz, σx, σy, σz)
    p = LVector(ωy = ωy, ωz = ωz, σx = σx, σy = σy, σz = σz)

    DiscreteDynamicalSystem(eom_verdes, u₀, p)
end

"""
    verdes(;u₀ = rand(3), ωy = 315, ωz = 80,
        σx = 0.0, σy = 0.0, σz = 0.0) → DiscreteDynamicalSystem

Intitialise a 3D system where the response X is a highly nonlinear combination
of Y and Z. The forcings Y and Z involve sines and cosines, and
have different periods, which controlled by `ωy` and `ωz`.

The equations of motion are

```math
\\begin{aligned}
x(t+1) &= \\dfrac{y(t)(18y(t) - 27y(t)^2 + 10)}{2} + z(t)(1-z(t)) + ηx \\
y(t+1) &= \\dfrac{(1 - \\dfrac{\\cos(2\\pi)}{\\omega y}t)}{2} + ηy \\
z(t+1) &= \\dfrac{(1 - \\dfrac{\\sin(2\\pi)}{\\omega z}t)}{2} + ηz
\\end{aligned}
```
where ηx, ηy, ηz is gaussian noise with mean 0 and standard deviation `σx`, `σy`
and `σz`.

# References
Verdes, P. F. "Assessing causality from multivariate time series." Physical
Review E 72.2 (2005): 026222.
"""
verdes(;u₀ = rand(3),
    ωy = 315, ωz = 80,
    σx = 0.01, σy = 0.01, σz = 0.01) =
    verdes(u₀, ωy, ωz, σx, σy, σz)
