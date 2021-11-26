using LabelledArrays

export logistic3

"""
    eom_logistic3(u, p, t)

Equations of motion for a system consisting of three coupled logistic map
representing the response of two independent dynamical variables to the
forcing from a common driver. The dynamical influence goes in the directions
Z → X and Z → Y.

## Equations of motion

The equations of motion are

```math
\\begin{aligned}
x(t+1) = (x(t)(r - r_1 x(t) - z(t) + σ_x η_x)) \\mod 1 \\\\
y(t+1) = (y(t)(r - r_2 y(t) - z(t) + σ_y η_y)) \\mod 1 \\\\
z(t+1) = (z(t)(r - r_3 z(t) + σ_z η_z)) \\mod 1
\\end{aligned}
```

Dynamical noise may be added to each of the dynamical variables by tuning the
parameters `σz`, `σx` and `σz`. Default values for the parameters
`r₁`, `r₂` and `r₃` are set such that the system exhibits chaotic behaviour,
with `r₁ = r₂ = r₃ = 4`.

## References

1. Runge, Jakob. Causal network reconstruction from time series: From theoretical 
    assumptions to practical estimation, Chaos 28, 075310 (2018); 
    doi: 10.1063/1.5025050
"""
function eom_logistic3(u, p, t)
    r₁, r₂, r₃, σx, σy, σz = (p...,)
    x, y, z = (u...,)

    # Independent dynamical noise for each variable.
    ηx = rand()
    ηy = rand()
    ηz = rand()

    dx = (x*(r₁ - r₁*x - z + σx*ηx)) % 1
    dy = (y*(r₂ - r₂*y - z + σy*ηy)) % 1
    dz = (z*(r₃ - r₃*z + σz*ηz)) % 1
    return SVector{3}(dx, dy, dz)
end

function logistic3(u₀, r₁, r₂, r₃, σx, σy, σz)
    p = @LArray [r₁, r₂, r₃, σx, σy, σz] (:r₁, :r₂, :r₃, :σx, :σy, :σz)
    DiscreteDynamicalSystem(eom_logistic3, u₀, p)
end

"""
    logistic3(;u₀ = rand(3), r = 4,
        σx = 0.05, σy = 0.05, σz = 0.05) → DiscreteDynamicalSystem

Initialise a dynamical system consisting of three coupled logistic map
representing the response of two independent dynamical variables to the
forcing from a common driver. The dynamical influence goes in the directions
``Z \\to X`` and ``Z \\to Y``.

## Equations of motion

The equations of motion are

```math
\\begin{aligned}
x(t+1) = (x(t)(r - r_1 x(t) - z(t) + σ_x η_x)) \\mod 1 \\\\
y(t+1) = (y(t)(r - r_2 y(t) - z(t) + σ_y η_y)) \\mod 1 \\\\
z(t+1) = (z(t)(r - r_3 z(t) + σ_z η_z)) \\mod 1
\\end{aligned}
```

Dynamical noise may be added to each of the dynamical variables by tuning the
parameters `σz`, `σx` and `σz`. Default values for the parameters
`r₁`, `r₂` and `r₃` are set such that the system exhibits chaotic behaviour,
with `r₁ = r₂ = r₃ = 4`.

## References

1. Runge, Jakob. Causal network reconstruction from time series: From theoretical 
    assumptions to practical estimation, Chaos 28, 075310 (2018); 
    doi: 10.1063/1.5025050
"""
logistic3(;u₀ = rand(3), r₁ = 4, r₂ = 4, r₃ = 4,
    σx = 0.05, σy = 0.05, σz = 0.05) = logistic3(u₀, r₁, r₂, r₃, σx, σy, σz)
