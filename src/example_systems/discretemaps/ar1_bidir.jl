using LabelledArrays

export ar1_bidir

"""
    eom_ar1_bidir(x, p, n) → SVector{2}

Equations of motion for a system consisting of two mutually 
coupled first order autoregressive processes. 

## Equations of motion 

```math
\\begin{aligned}
x(t+1) &= a_{1}x + c_{yx}y + \\epsilon_{x} \\\\
y(t+1) &= b_{1}y + c_{xy}x + \\epsilon_{y}
\\end{aligned}
```
    
where ``\\epsilon_{x}`` and ``\\epsilon_{y}`` are drawn independently 
at each time step from normal distributions with zero mean and standard 
deviations `σx` and `σy`.
"""
function eom_ar1_bidir(x, p, n)
    a₁, b₁, c_xy, c_yx, σx, σy = (p...,)
    x, y = (x...,)

    ϵx = rand(Normal(0, σx))
    ϵy = rand(Normal(0, σy))

    dx = a₁*x + c_yx*y + ϵx
    dy = b₁*y + c_xy*x + ϵy
    return SVector{2}(dx, dy)
end

function ar1_bidir(u₀,a₁, b₁, c_xy, c_yx, σx, σy)
    p = @LArray [a₁, b₁, c_xy, c_yx, σx, σy] (:a₁, :b₁, :c_xy, :c_yx, :σx, :σy)
    logistic_system = DiscreteDynamicalSystem(eom_ar1_bidir, u₀, p)

    return logistic_system
end

"""
    ar1_bidir(;u₀ = rand(2), a₁ = 0.5, b₁ = 0.7, c_xy = 0, c_yx = 0.2, 
        σx = 0.3, σy = 0.3) → DiscreteDynamicalSystem

A system consisting of two mutually coupled first order autoregressive processes. 

## Equations of motion 

```math
\\begin{aligned}
x(t+1) &= a_{1}x + c_{yx}y + \\epsilon_{x} \\\\
y(t+1) &= b_{1}y + c_{xy}x + \\epsilon_{y}
\\end{aligned}
```

where ``\\epsilon_{x}`` and ``\\epsilon_{y}`` are drawn independently 
at each time step from normal distributions with zero mean and standard 
deviations `σx` and `σy`.
"""
ar1_bidir(;a₁ = 0.5, b₁ = 0.7, u₀ = rand(2), c_xy = 0, c_yx = 0.2, σx = 0.3, σy = 0.3) =
    ar1_bidir(u₀, a₁, b₁, c_xy, c_yx, σx, σy)
