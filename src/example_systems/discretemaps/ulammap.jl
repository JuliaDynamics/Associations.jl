using LabelledArrays 

export ulam

function eom_ulam(dx, x, p, t)
    ε = p[:ε]
    f = x -> 2 - x^2
    dx[1] = f(ε*x[length(dx)] + (1-ε)*x[1])
    for i in 2:length(dx)
        dx[i] = f(ε*x[i-1] + (1-ε)*x[i])
    end
end

"""
    ulam(D::Int = 10; u₀ = rand(D), ε::Real = 0.10) → DiscreteDynamicalSystem

A lattice of `D` unidirectionally coupled ulam maps[^Schreiber2000] defined as 

```math
x^{m}_{t+1} = f(\\epsilon x^{m-1}_{t} + (1 - \\epsilon) x_{t}^{m}),
```

where ``m = 1, 2, \\ldots, D`` and ``f(x) = 2 - x^2``. In this system, information transfer 
happens only in the direction of increasing ``m``.

[^Schreiber2000]: Schreiber, Thomas. "Measuring information transfer." Physical review letters 85.2 (2000): 461.
"""
function ulam(D::Int = 10; u₀ = rand(D), ε::Real = 0.10)

    p = LVector(ε = ε)

    DiscreteDynamicalSystem(eom_ulam, u₀, p)
end