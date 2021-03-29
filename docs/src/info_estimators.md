# [Entropy/probability estimators](@id estimators)

Information theoretic causality measures in this package are calculated using entropy estimation. To do so, it uses estimators and [`genentropy`](@ref) from the [Entropies.jl](https://github.com/JuliaDynamics/Entropies.jl) package. However, additional estimators are also available for some of the higher-level methods.

An exhaustive list of estimators used by either [`transferentropy`](@ref), [`predictive_asymmetry`](@ref), [`mutualinfo`](@ref) or [`genentropy`](@ref) is given below.

*Note: Not all estimators are valid for all methods. Check individual method docstrings for more information.*

## Binning based

### Visitation frequency

```@docs
VisitationFrequency
RectangularBinning
```

### Transfer operator

```@docs
TransferOperator
```

## Kernel density based

```@docs
NaiveKernel
TreeDistance
DirectDistance
```

## Nearest neighbor based

```@docs
KozachenkoLeonenko
Kraskov
Kraskov1
Kraskov2
```

## Permutation based

```@docs
SymbolicPermutation
```

## Hilbert

```@docs
Hilbert
Amplitude
Phase
```
