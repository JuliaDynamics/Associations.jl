# [Estimators](@id estimators)

Information theoretic causality measures in this package are calculated using entropy estimation. To do so, it uses estimators and [`genentropy`](@ref) from the [Entropies.jl](https://github.com/JuliaDynamics/Entropies.jl) package. However, additional estimators are also available for some of the higher-level methods.

The following method-estimator combinations are possible.

| Estimator                       | [`genentropy`](@ref) | [`mutualinfo`](@ref) | [`transferentropy`](@ref) | [`predictive_asymmetry`](@ref) |
|:------------------------------- | :-: |:-: |:-: |:-: |
| [`VisitationFrequency`](@ref)                        | ✓ | ✓ | ✓ | ✓ |
| [`TransferOperator`](@ref)                           | ✓ | ✓ | ✓ | ✓ |
| [`NaiveKernel`](@ref)                                | ✓ | ✓ | ✓ | ✓ |
| [`SymbolicPermutation`](@ref)                        | ✓ | ✓ | ✓ | ✓ |
| [`SymbolicAmplitudeAwarePermutation`](@ref symbolic) | ✓ |   |   |   |
| [`SymbolicWeightedPermutation`](@ref symbolic)       | ✓ |   |   |   |
| [`KozachenkoLeonenko`](@ref)                         | ✓ | ✓ | ✓ | ✓ |
| [`Kraskov`](@ref)                                    | ✓ | ✓ | ✓ | ✓ |
| [`Kraskov1`](@ref)                                   | ✓ | ✓ |   |   |
| [`Kraskov2`](@ref)                                   | ✓ | ✓ |   |   |
| [`Hilbert`](@ref)                                    | ✓ |   | ✓ | ✓ |
| [`TimeScaleMODWT`](@ref)                             | ✓ |   |   |   |

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

## [Permutation based](@id symbolic)

```@docs
SymbolicPermutation
```

## Wavelet based

```@docs
TimeScaleMODWT
```

## Hilbert

```@docs
Hilbert
Amplitude
Phase
```
