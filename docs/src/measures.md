# [Association measures](@id association_measure)

## Overview

| Type                 | Measure                               | Pairwise | Conditional | Function version               |
| -------------------- | ------------------------------------- | :------: | :---------: | ------------------------------ |
| Correlation          | [`PearsonCorrelation`](@ref)          |    ✓    |     ✖      | [`pearson_correlation`](@ref)  |
| Correlation          | [`DistanceCorrelation`](@ref)         |    ✓    |     ✖      | [`distance_correlation`](@ref) |
| Closeness            | [`SMeasure`](@ref)                    |    ✓    |     ✖      | [`s_measure`](@ref)            |
| Closeness            | [`HMeasure`](@ref)                    |    ✓    |     ✖      | [`h_measure`](@ref)            |
| Closeness            | [`MMeasure`](@ref)                    |    ✓    |     ✖      | [`m_measure`](@ref)            |
| Closeness (ranks)    | [`LMeasure`](@ref)                    |    ✓    |     ✖      | [`l_measure`](@ref)            |
| Closeness            | [`JointDistanceDistribution`](@ref)   |    ✓    |     ✖      | [`jdd`](@ref)                  |
| Cross-mapping        | [`PairwiseAsymmetricInference`](@ref) |    ✓    |     ✖      | [`crossmap`](@ref)             |
| Cross-mapping        | [`ConvergentCrossMapping`](@ref)      |    ✓    |     ✖      | [`crossmap`](@ref)             |
| Shared information   | [`MIShannon`](@ref)                   |    ✓    |     ✖      | [`mutualinfo`](@ref)           |
| Shared information   | [`MIRenyiJizba`](@ref)                |    ✓    |     ✖      | [`mutualinfo`](@ref)           |
| Shared information   | [`MIRenyiSarbu`](@ref)                |    ✓    |     ✖      | [`mutualinfo`](@ref)           |
| Shared information   | [`MITsallisFuruichi`](@ref)           |    ✓    |     ✖      | [`mutualinfo`](@ref)           |
| Shared information   | [`PartialCorrelation`](@ref)          |    ✖    |     ✓      | [`partial_correlation`](@ref)  |
| Shared information   | [`CMIShannon`](@ref)                  |    ✖    |     ✓      | [`condmutualinfo`](@ref)       |
| Shared information   | [`CMIRenyiSarbu`](@ref)               |    ✖    |     ✓      | [`condmutualinfo`](@ref)       |
| Shared information   | [`CMIRenyiJizba`](@ref)               |    ✖    |     ✓      | [`condmutualinfo`](@ref)       |
| Information transfer | [`TEShannon`](@ref)                   |    ✓    |     ✓      | [`transferentropy`](@ref)      |
| Information transfer | [`TERenyiJizba`](@ref)                |    ✓    |     ✓      | [`transferentropy`](@ref)      |

## Correlation measures

### Pearson correlation

```@docs
PearsonCorrelation
pearson_correlation
```

### Partial correlation

```@docs
PartialCorrelation
partial_correlation
```

### Distance correlation

```@docs
DistanceCorrelation
distance_correlation
```

## Closeness measures

### Joint distance distribution

```@docs
JointDistanceDistribution
jdd
```

### S-measure

```@docs
SMeasure
s_measure
```

### H-measure

```@docs
HMeasure
h_measure
```

### M-measure

```@docs
MMeasure
m_measure
```

### L-measure

```@docs
LMeasure
l_measure
```

## Cross-map measures

See also the [cross mapping API](@ref crossmap_api) for estimators.

### Convergent cross mapping

```@docs
ConvergentCrossMapping
```

### Pairwise asymmetric inference

```@docs
PairwiseAsymmetricInference
```

## [Information measures](@id information_measures)

Association measures that are information-based are listed here. Available estimators
are listed in the [information API](@ref information_api).

### Mutual information (Shannon)

```@docs
MIShannon
```

### Mutual information (Tsallis, Furuichi)

```@docs
MITsallisFuruichi
```

### Mutual information (Tsallis, Martin)

```@docs
MITsallisMartin
```

### Mutual information (Rényi, Sarbu)

```@docs
MIRenyiSarbu
```

### Mutual information (Rényi, Jizba)

```@docs
MIRenyiJizba
```

### Conditional mutual information (Shannon)

```@docs
CMIShannon
```

### Conditional mutual information (Rényi, Jizba)

```@docs
CMIRenyiJizba
```

### Conditional mutual information (Rényi, Poczos)

```@docs
CMIRenyiPoczos
```

### Transfer entropy (Shannon)

```@docs
TEShannon
```

### Transfer entropy (Rényi, Jizba)

```@docs
TERenyiJizba
```
