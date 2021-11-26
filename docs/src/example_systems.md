# Example systems

## [Continuous](@id continuous_systems)

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
lorenz_lorenz_lorenz_transitive
```

### Two bidirectionally coupled 3D Rössler systems

```@docs
rossler_rossler_bidir
```

### Two bidirectionally coupled 3D Rössler systems forced by another 3D Rössler system

```@docs
rossler_rossler_rossler_bidir_forced
```

### Unidirectonal forcing from a 3D Rössler system to a 3D Lorenz system

```@docs
rossler_lorenz
```

### N-scroll Chua attractors

```@docs
chuacircuit_nscroll_sine
```

## [Discrete](@id discrete_systems)

### [Ulam map](@id ulam)

```@docs
ulam
```

### [Autoregressive order one 2D system](@id system_ar1)

```@docs
ar1_unidir
```

### [Nonlinear 3D system with nonlinear coupling](@id system_nonlinear3d)

```@docs
nonlinear3d
```

### [Unidirectionally coupled 2D logistic maps](@id system_logistic2_unidir)

```@docs
logistic2_unidir
```

### [Bidirectionally coupled 2D logistic maps](@id system_logistic2_bidir)

```@docs
logistic2_bidir
```

### Forcing of two independent logistic maps from common logistic map driver

```@docs
logistic3
```

### Unidirectional, transitive chain of logistic maps

```@docs
logistic4
```

### [Two unidirectionally coupled Henon maps](@id henon2d)

```@docs
henon2
```

### Strange, nonchaotic attractors

```@docs
anishchenko1
```

### [Bidirectional Ikeda map](@id ikeda2d_bidir)

```@docs
ikeda
```

### Noise

```@docs
noise_uu
```

```@docs
noise_ug
```

```@docs
noise_brownian
```
