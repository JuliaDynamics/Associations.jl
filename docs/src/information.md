# Information API

```@docs
CausalityTools.information(::MultivariateInformationMeasureEstimator)
```

## Definitions

We implement a range of bivariate and multivariate information measures (i.e. measures
that are functionals of probability mass functions or probability densities). They are 
listed below. For estimating a measure, use [`information`](@ref) with a compatible
estimator among the estimators listed below.

### [Conditional entropies](@id conditional_entropies)

```@docs
ConditionalEntropy
ConditionalEntropyShannon
ConditionalEntropyTsallisFuruichi
ConditionalEntropyTsallisAbe
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
MutualInformation
MIShannon
MITsallisFuruichi
MITsallisMartin
MIRenyiJizba
MIRenyiSarbu
```

### Conditional mutual informations

```@docs
ConditionalMutualInformation
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
MutualInformationEstimator
KraskovStögbauerGrassberger1
KraskovStögbauerGrassberger2
GaoKannanOhViswanath
GaoOhViswanath
GaussianMI
```


### Conditional mutual information estimators

```@docs
ConditionalMutualInformationEstimator
GaussianCMI
FPVP
MesnerShalizi
Rahimzamani
PoczosSchneiderCMI
```


## Convenience functions

For commonly used names, we provide convenience functions. These are just simple 
wrappers around [`information`](@ref).

```@docs
joint_entropy
conditional_entropy
mutualinfo
condmutualinfo
```

## Single-variable information API

```@docs

```

