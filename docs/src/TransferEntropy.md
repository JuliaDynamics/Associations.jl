# [Transfer entropy](@ref transferentropy)

The [`transferentropy`](@ref) and [`mutualinfo`](@ref) functions uses estimators from the
[Entropies.jl](https://github.com/JuliaDynamics/Entropies.jl) to compute transfer entropy
and mutual information, respectively.

```@docs
transferentropy
```

## [Estimators](@id estimators)

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
