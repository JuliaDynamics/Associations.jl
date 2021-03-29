# Example systems

## [Continuous coupled dynamical systems](@id continuous_systems)

### Mediated link

```@docs
ExampleSystems.mediated_link(;u₀ = rand(9), ωx = 1, ωy = 1.015, ωz = 0.985,
    k = 0.15, l = 0.2, m = 10.0, c = 0.06)
```

### Two bidirectionally coupled 3D Lorenz systems

```@docs
ExampleSystems.lorenz_lorenz_bidir
```

### Two bidirectionally coupled 3D Lorenz systems forced by another 3D Lorenz system

```@docs
lorenz_lorenz_lorenz_bidir_forced
```

### Three transitively connected 3D Lorenz systems

```@docs
lorenz_lorenz_lorenz_transitive(;uᵢ=rand(9),
        σ₁ = 10.0, σ₂ = 10.0, σ₃ = 10.0,
        ρ₁ = 28.0, ρ₂ = 28.0, ρ₃ = 28.0,
        β₁ = 8/3,  β₂ = 8/3,  β₃ = 8.3,
        c₁₂ = 1.0, c₂₃ = 1.0)
```

### Two bidirectionally coupled 3D Rössler systems

```@docs
rossler_rossler_bidir(; u0 = rand(6),
        ω₁ = 1.015, ω₂ = 0.985,
        c_xy = 0.1, c_yx = 0.1,
        a₁ = 0.15, a₂ = 0.2, a₃ = 10,
        b₁ = 0.15, b₂ = 0.2, b₃ = 10)
```

### Two bidirectionally coupled 3D Rössler systems forced by another 3D Rössler system

```@docs
rossler_rossler_rossler_bidir_forced(; u0 = rand(9),
        ω₁ = 1.015, ω₂ = 0.985, ω₃ = 0.95,
        c_xy = 0.1, c_yx = 0.1,
        c_zx = 0.05, c_zy = 0.05,
        a₁ = 0.15, a₂ = 0.2, a₃ = 10,
        b₁ = 0.15, b₂ = 0.2, b₃ = 10,
        c₁ = 0.15, c₂ = 0.2, c₃ = 10)
```

### Unidirectonal forcing from a 3D Rössler system to a 3D Lorenz system

```@docs
rossler_lorenz(;u₀ = rand(6), a₁ = -6, a₂ = 6, a₃ = 2.0,
        b₁ = 10, b₂ = 28, b₃ = 8/3, c_xy = 1)
```

### N-scroll Chua attractors

```@docs
chuacircuit_nscroll_sine(;u₀ = [0.0, 0.0, 0.28695],
        α = 10.814, β = 14, γ = 0, a = 1.3, b = 0.11, c = 2,
        σx = 0.0, σy = 0.0, σz = 0.0)
```

## [Discrete coupled dynamical systems](@id discrete_systems)

### [Ulam map](@id ulam)

```@docs
ulam
```

### [Autoregressive order one 2D system](@id system_ar1)

```@docs
ar1_unidir(;uᵢ = rand(2), a₁ = 0.90693, b₁ = 0.40693, c_xy = 0.5, σ = 0.40662)
```

### [Nonlinear 3D system with nonlinear coupling](@id system_nonlinear3d)

```@docs
nonlinear3d
```

### [Unidirectionally coupled 2D logistic maps](@id system_logistic2_unidir)

```@docs
logistic2_unidir(;u₀ = rand(2), c_xy = 0.1, r₁ = 3.78, r₂ = 3.66, σ = 0.05)
```

### [Bidirectionally coupled 2D logistic maps](@id system_logistic2_bidir)

```@docs
logistic2_bidir(;u₀ = rand(2), c_xy = 0.1, c_yx = 0.1,
    r₁ = 3.78, r₂ = 3.66, σ_xy = 0.05, σ_yx = 0.05)
```

### Forcing of two independent logistic maps from common logistic map driver

```@docs
logistic3(;u₀ = rand(3), r₁ = 4, r₂ = 4, r₃ = 4, σx = 0.05, σy = 0.05, σz = 0.05)
```

### Unidirectional, transitive chain of logistic maps

```@docs
logistic4(;u₀ = rand(4),
            r₁ = 3.9, r₂ = 3.6, r₃ = 3.6, r₄ = 3.8,
            c₁₂ = 0.4, c₂₃ = 0.4, c₃₄ = 0.35)
```

### [Two unidirectionally coupled Henon maps](@id henon2d)

```@docs
henon2(;u₀ = rand(4), c_xy = 2.0)
```

### Strange, nonchaotic attractors

```@docs
anishchenko1(;u₀ = rand(2), α =3.277, s=0.1, ω=0.5*(sqrt(5)-1))
```

### [Bidirectional Ikeda map](@id ikeda2d_bidir)

```@docs
ikeda
```

## Noise

```@docs
noise_uu
```

```@docs
noise_ug
```

```@docs
noise_brownian
```
