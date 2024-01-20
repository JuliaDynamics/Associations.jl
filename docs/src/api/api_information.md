# [Information API](@ref information_api)

The information API define multivariate information measures, which we define as 
any functional of probability mass functions (pmf) or probability densities that operate
on two or more input datasets.

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

## Core

### Discretization

```@docs
Discretization
codify
```

#### Column-wise discretization

```@docs
CodifyVariables
```

#### Point-wise discretization

```@docs
CodifyPoints
Encoding
encoding
GaussianCDFEncoding
OrdinalPatternEncoding
RelativeMeanEncoding
RelativeFirstDifferenceEncoding
UniqueElementsEncoding
CombinationEncoding
RectangularBinEncoding
encode
decode
```

## Counting and probabilities

For counting and probabilities, CausalityTools.jl extends the single-variable machinery
in ComplexityMeasures.jl to multiple variables.

```@docs
CausalityTools.Counts
CausalityTools.counts
CausalityTools.Probabilities
CausalityTools.probabilities
marginal
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
CMIDecomposition
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

### Transfer entropy estimators

```@docs
TransferEntropyEstimator
Zhu1
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

## Single-variable information API (from ComplexityMeasures.jl)

Below we list some relevant types from
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl) that 
are used for the [`EntropyDecomposition`](@ref) estimator.

### Entropies

```@docs
Shannon
Renyi
Tsallis
```

### Discrete information estimators

```@docs
DiscreteInfoEstimator
PlugIn
MillerMadow
Schuermann
GeneralizedSchuermann
Jackknife
HorvitzThompson
ChaoShen
```

```@docs
OutcomeSpace
```

#### Binning

```@docs
ValueBinning
RectangularBinning
FixedRectangularBinning
```

#### Ordinal patterns

```@docs
OrdinalPatterns
```

#### Dispersion

```@docs
Dispersion
```

#### Unique elements

The [`UniqueElements`](@ref) outcome space is useful for categorical data.

```@docs
UniqueElements
```

### Differential information estimators

```@docs
DifferentialInfoEstimator
Kraskov
KozachenkoLeonenko
Zhu
ZhuSingh
Gao
Goria
Lord
LeonenkoProzantoSavani
Vasicek
AlizadehArghami
Ebrahimi
Correa
```
