
# [Discrete coupled dynamical systems](@id discrete_systems)

## Autoregressive order one 2D system

```@docs
ar1_unidir(;uᵢ = rand(2), a₁ = 0.90693, b₁ = 0.40693, c_xy = 0.5, σ = 0.40662)
```

## Nonlinear 3D system with nonlinear coupling

```@docs
nonlinear3d
```

## Unidirectionally coupled 2D logistic maps

```@docs
logistic2_unidir(;u₀ = rand(2), c_xy = 0.1, r₁ = 3.78, r₂ = 3.66, σ = 0.05)
```

## Bidirectionally coupled 2D logistic maps

```@docs
logistic2_bidir(;u₀ = rand(2), c_xy = 0.1, c_yx = 0.1,
    r₁ = 3.78, r₂ = 3.66, σ_xy = 0.05, σ_yx = 0.05)
```

## Forcing of two independent logistic maps from common logistic map driver

```@docs
logistic3(;u₀ = rand(3), r₁ = 4, r₂ = 4, r₃ = 4, σx = 0.05, σy = 0.05, σz = 0.05)
```

## Unidirectional, transitive chain of logistic maps

```@docs
logistic4(;u₀ = rand(4),
            r₁ = 3.9, r₂ = 3.6, r₃ = 3.6, r₄ = 3.8,
            c₁₂ = 0.4, c₂₃ = 0.4, c₃₄ = 0.35)
```

## [Two unidirectionally coupled Henon maps](@id henon2d)

```@docs
henon2(;u₀ = rand(4), c_xy = 2.0)
```

## Strange, nonchaotic attractors

```@docs 
anishchenko1(;u₀ = rand(2), α =3.277, s=0.1, ω=0.5*(sqrt(5)-1))
```
