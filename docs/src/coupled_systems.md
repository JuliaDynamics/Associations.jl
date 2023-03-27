# Predefined coupled systems

## Systems definition API

The systems definition API is defined by

- [`SystemDefinition`](@ref), [`DiscreteDefinition`](@ref), and [`ContinuousDefinition`](@ref).
- [`system`](@ref)

```@docs
SystemDefinition
DiscreteDefinition
ContinuousDefinition
system
CausalityTools.SimpleDiGraph
CausalityTools.SimpleWeightedDiGraph
```

## Discrete systems

```@docs
Anishchenko
AR1Unidir
AR1Bidir
ChaoticMaps3
Henon2
Henon3
Ikeda2
ChaoticNoisyLinear2
Logistic2Unidir
Logistic2Bidir
Logistic3CommonDriver
Logistic4Chain
Nonlinear3
Peguin2
UlamLattice
Var1
Verdes3
```

## Continuous systems

```@docs
ChuaCircuitsBidir6
ChuaScrollSine3
HindmarshRose3
LorenzBidir6
LorenzForced9
MediatedLink9
Repressilator6
RosslerBidir6
RosslerForced9
RosslerLorenzUnidir6
Thomas3
```
