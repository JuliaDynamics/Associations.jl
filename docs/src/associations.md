# Association measures

The most basic components of CausalityTools.jl are a collection of statistics that in some manner quantify the "association" between input datasets. Precisely what is meant by "association" depends on the measure, and precisely what is meant  by "quantify" depends on the *estimator* of that measure. 

Many association statistics exists. We have divided these statistics into (currently) three separate APIs:

- The [information API](@ref information_api).
- The [cross mapping API](@ref cross_mapping_api).
- The [correlation API](@ref correlation_api).

## [Information API](@ref information_api)

The information API defines bivariate and multivariate information measures, which we here define as any functional of a multidimensional probability mass functions (PMFs) or multidimensional probability density.

### Overview

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

### Discretization

A fundamental operation when computing multivariate information measures from data is *discretization*. There are many ways of discretizing multiple input datasets. We offer two main ways of doing so.

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

### Counting and probabilities

For counting and probabilities, CausalityTools.jl extends the single-variable machinery
in ComplexityMeasures.jl to multiple variables.

```@docs
CausalityTools.Counts
CausalityTools.counts
CausalityTools.Probabilities
CausalityTools.probabilities
marginal
```


### [Definitions](@id definitions)


#### [Conditional entropies](@id conditional_entropies)

```@docs
ConditionalEntropy
ConditionalEntropyShannon
ConditionalEntropyTsallisFuruichi
ConditionalEntropyTsallisAbe
```

#### [Divergences and distances](@id divergences_and_distances)

```@docs
HellingerDistance
KLDivergence
RenyiDivergence
VariationDistance
```

#### [Joint entropies](@id joint_entropies)

```@docs
JointEntropy
JointEntropyShannon
JointEntropyTsallis
JointEntropyRenyi
```

#### Mutual informations

```@docs
MutualInformation
MIShannon
MITsallisFuruichi
MITsallisMartin
MIRenyiJizba
MIRenyiSarbu
```

#### Conditional mutual informations

```@docs
ConditionalMutualInformation
CMIShannon
CMIRenyiSarbu
CMIRenyiJizba
CMIRenyiPoczos
CMITsallis
```

#### Partial mutual information

```@docs
PartialMutualInformation
```

#### Transfer entropy

```@docs
TransferEntropy
TEShannon
TERenyiJizba
```

### Estimators

#### Generic estimators

```@docs
JointProbabilities
EntropyDecomposition
MIDecomposition
CMIDecomposition
```

#### Mutual information estimators

```@docs
MutualInformationEstimator
KraskovStögbauerGrassberger1
KraskovStögbauerGrassberger2
GaoKannanOhViswanath
GaoOhViswanath
GaussianMI
```


#### Conditional mutual information estimators

```@docs
ConditionalMutualInformationEstimator
GaussianCMI
FPVP
MesnerShalizi
Rahimzamani
PoczosSchneiderCMI
```

#### Transfer entropy estimators

```@docs
TransferEntropyEstimator
Zhu1
```

### [Convenience functions](@ref convenience_info)

For commonly used names, we provide convenience functions. These are just simple 
wrappers around [`information`](@ref).

```@docs
joint_entropy
conditional_entropy
mutualinfo
condmutualinfo
```

### Single-variable information API (from ComplexityMeasures.jl)

Below we list some relevant types from
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl) that 
are used for the [`EntropyDecomposition`](@ref) estimator.

#### Entropies

```@docs
Shannon
Renyi
Tsallis
Kaniadakis
```

#### Discrete information estimators

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

#### Differential information estimators

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

## [Correlation API](@id correlation_api)

This page lists all available [`CorrelationMeasure`](@ref)s, as 
well as their convenience functions. The [examples](@ref correlation_examples)
is also useful.

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

## [Cross mapping API](@id cross_mapping_api)

The cross mapping API define different ways of quantifying association based on the 
concept of "cross mapping", which has appeared in many contexts in the literature,
and gained huge popularity with  [Sugihara2012](@citet)'s on *convergent cross mapping*.

Since their paper, several cross mapping methods and frameworks have emerged in the
literature. In CausalityTools.jl, we provide a unified interface for using these cross
mapping methods.

To estimate a cross map measure, you simply input a [`CrossmapMeasure`](@ref) instance
as the first argument to a [`CrossmapEstimator`](@ref), which is then fed to 
the [`crossmap`](@ref) or [`predict`](@ref) functions. 

The cross mapping API consists of the following functions.

- [`predict`](@ref)
- [`crossmap`](@ref)

These functions can dispatch on a [`CrossmapMeasure`](@ref), and we currently implement

- [`ConvergentCrossMapping`](@ref).
- [`PairwiseAsymmetricEmbedding`](@ref).

```@docs
crossmap
predict
```

### Measures

```@docs
CrossmapMeasure
ConvergentCrossMapping
PairwiseAsymmetricInference
```

### Estimators

```@docs
CrossmapEstimator
RandomVectors
RandomSegment
ExpandingSegment
```
