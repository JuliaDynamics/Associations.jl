# Information API

```@docs
CausalityTools.information
```

## Definitions

We implement a range of bivariate and multivariate information measures (i.e. measures
that are functionals of probability mass functions or probability densities). They are 
listed below. For estimating a measure, use [`information`](@ref) with a compatible
estimator among the estimators listed below.

### [Conditional entropies](@id conditional_entropies)

```@docs
CEShannon
CETsallisAbe
CETsallisFuruchi
```

### [Divergences and distances](@id divergences_and_distances)

```@docs
HellingerDistance
KLDivergence
RenyiDivergence
VariationDistance
```

### [Joint entropies](@id joint_entropies)

```@docs
JointEntropyShannon
JointEntropyTsallis
JointEntropyRenyi
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

### Other measures

```@docs
PartialMutualInformation
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

## Single-variable information API

```@docs
Goria
```