# Information API

```@docs
information
```

## Definitions

We implement a range of bivariate and multivariate information measures (i.e. measures
that are functionals of probability mass functions or probability densities). They are 
listed below. For estimating a measure, use [`information`](@ref) with a compatible
estimator among the estimators listed below.

### Conditional entropies

```@docs
CEShannon
CETsallisAbe
CETsallisFuruchi
```

### Divergences and distances

```@docs
HellingerDistance
KLDivergence
RenyiDivergence
VariationDistance
```

### Mutual informations

```@docs
MIShannon
MITsallisFuruichi
MITsallisMartin
MIRenyiJizba
MIRenyiSarbu
```

### Conditional mutual informations

```@docs
CMIShannon
CMIRenyiSarbu
CMIRenyiJizba
CMIRenyiPoczos
CMITsallis
```

## Estimators

### Generic estimators
```@docs
JointProbabilities
EntropyDecomposition
MIDecomposition
```

### Mutual information estimators

```@docs
KSG1
KSG2
GaoKannanOhViswanath
GaoOhViswanath
GaussianMI
```

