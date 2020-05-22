
# [Continuous coupled dynamical systems](@id continuous_systems)

## Mediated link

```@docs
mediated_link(;u₀ = rand(9), ωx = 1, ωy = 1.015, ωz = 0.985,
    k = 0.15, l = 0.2, m = 10.0, c = 0.06)
```

## Two bidirectionally coupled 3D Lorenz systems

```@docs
lorenz_lorenz_bidir(; u0 = rand(6),
        c_xy = 0.2, c_yx = 0.2,
        a₁ = 10, a₂ = 28, a₃ = 8/3,
        b₁ = 10, b₂ = 28, b₃ = 9/3)
```

## Two bidirectionally coupled 3D Lorenz systems forced by another 3D Lorenz system

```@docs
lorenz_lorenz_lorenz_bidir_forced(; u0 = rand(9),
        c_xy = 0.1, c_yx = 0.1,
        c_zx = 0.05, c_zy = 0.05,
        a₁ = 10, a₂ = 28, a₃ = 8/3,
        b₁ = 10, b₂ = 28, b₃ = 8/3,
        c₁ = 10, c₂ = 28, c₃ = 8/3)
```

## Three transitively connected 3D Lorenz systems

```@docs
lorenz_lorenz_lorenz_transitive(;uᵢ=rand(9),
        σ₁ = 10.0, σ₂ = 10.0, σ₃ = 10.0,
        ρ₁ = 28.0, ρ₂ = 28.0, ρ₃ = 28.0,
        β₁ = 8/3,  β₂ = 8/3,  β₃ = 8.3,
        c₁₂ = 1.0, c₂₃ = 1.0)
```

## Two bidirectionally coupled 3D Rössler systems

```@docs
rossler_rossler_bidir(; u0 = rand(6),
        ω₁ = 1.015, ω₂ = 0.985,
        c_xy = 0.1, c_yx = 0.1,
        a₁ = 0.15, a₂ = 0.2, a₃ = 10,
        b₁ = 0.15, b₂ = 0.2, b₃ = 10)
```

## Two bidirectionally coupled 3D Rössler systems forced by another 3D Rössler system

```@docs
rossler_rossler_rossler_bidir_forced(; u0 = rand(9),
        ω₁ = 1.015, ω₂ = 0.985, ω₃ = 0.95,
        c_xy = 0.1, c_yx = 0.1,
        c_zx = 0.05, c_zy = 0.05,
        a₁ = 0.15, a₂ = 0.2, a₃ = 10,
        b₁ = 0.15, b₂ = 0.2, b₃ = 10,
        c₁ = 0.15, c₂ = 0.2, c₃ = 10)
```

## Unidirectonal forcing from a 3D Rössler system to a 3D Lorenz system

```@docs
rossler_lorenz(;u₀ = rand(6), a₁ = -6, a₂ = 6, a₃ = 2.0,
        b₁ = 10, b₂ = 28, b₃ = 8/3, c_xy = 1)
```

## N-scroll Chua attractors

```@docs
chuacircuit_nscroll_sine(;u₀ = [0.0, 0.0, 0.28695],
        α = 10.814, β = 14, γ = 0, a = 1.3, b = 0.11, c = 2,
        σx = 0.0, σy = 0.0, σz = 0.0)
```
