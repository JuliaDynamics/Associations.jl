
function eom_latticeunidir(dx, x, p, t)
    f, ε = p[1], p[2]
    dx[1] = f(ε*x[length(dx)] + (1-ε)*x[1])
    for i in 2:length(dx)
        dx[i] = f(ε*x[i-1] + (1-ε)*x[i])
    end
end

"""
    latticeunidir(D::Int; uᵢ = rand(D), ϵ::Real = 0.10, f::Function = x -> 2 - x^2) → DiscreteDynamicalSystem

A lattice of `D` unidirectionally coupled maps[^Schreiber2000] defined as 

```math
x^{m}_{t+1} = f(\\epsilon x^{m-1}_{t} + (1 - \\epsilon) x_{t}^{m}),
```

where ``m = 1, 2, \\ldots, D``. In this system, information transfer happens only in the direction of increasing ``m``.

The default `f`, ``f(x) = 2 - x^2``, is the Ulam map used in Schreiber (2000). 

[^Schreiber2000]: Schreiber, Thomas. "Measuring information transfer." Physical review letters 85.2 (2000): 461.
"""
function latticeunidir(D::Int; uᵢ = rand(D), ϵ::Real = 0.10, f::Function = x -> 2 - x^2)
    p = (f, ϵ)
    DiscreteDynamicalSystem(eom_latticeunidir, uᵢ, p)
end