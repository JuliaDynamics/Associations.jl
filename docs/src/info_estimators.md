## [Entropy/probability estimators](@id estimators)

An exhaustive list of estimators used by either [`transferentropy`](@ref), [`predictive_asymmetry`](@ref), [`mutualinfo`](@ref) or [`genentropy`](@ref). 

*Note: Not all estimators are valid for all methods. Check individual method docstrings for more information.*

### Binning based

#### Visitation frequency

```@docs
VisitationFrequency
RectangularBinning
```

#### Transfer operator

```@docs
TransferOperator
invariantmeasure
InvariantMeasure
transfermatrix
```

### Kernel density based

```@docs
NaiveKernel
TreeDistance
DirectDistance
```

### Nearest neighbor based

```@docs
KozachenkoLeonenko
Kraskov
Kraskov1
Kraskov2
```

### Permutation based

```@docs
SymbolicPermutation
```

### Hilbert

```@docs
Hilbert
Amplitude
Phase
```
