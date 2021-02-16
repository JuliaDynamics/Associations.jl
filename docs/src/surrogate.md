# Surrogate data

Surrogate time series can be used for hypothesis testing in time series causality analyses. 

The [TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/Entropies.jl) package provides methods for generating surrogate time series. Relevant methods are re-exported here for convenience. For more details, see its package documentation.

## Generation of surrogate data

TimeseriesSurrogates.jl exports two main functions. Both of them dispatch on the chosen method, a subtype of `Surrogate`.

```@docs
surrogate
surrogenerator
```

## Surrogate methods

```@docs
RandomShuffle
BlockShuffle
CycleShuffle
CircShift
RandomFourier
TFTS
AAFT
TAAFT
IAAFT
AutoRegressive
PseudoPeriodic
WLS
ShuffleDimensions
```

### Utils

```@docs
noiseradius
```
