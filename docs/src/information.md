# Information API

## Overview

CausalityTools.jl implements a range of bivariate and multivariate information measures,
which are listed under [definition](@ref definitions) below. Each of these measures
can be estimated using one or several [estimators](@ref) with [`information`](@ref).
Please see the [tutorial](@ref info_tutorial) for examples.

The API consists of the following types and methods:
- [`MultivariateInformationMeasure`](@ref)
- [`MultivariateInformationMeasureEstimator`](@ref)
- [`information`](@ref)

```@docs
MultivariateInformationMeasure
MultivariateInformationMeasureEstimator
CausalityTools.information(::MultivariateInformationMeasureEstimator)
```

## [Definitions](@id definitions)


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
JointEntropy
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

### Partial mutual information

```@docs
PartialMutualInformation
```

### Transfer entropy

```@docs
TransferEntropy
TEShannon
TERenyiJizba
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


## [Convenience functions](@ref convenience_info)

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

