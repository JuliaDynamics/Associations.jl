# Association measures

The most basic components of CausalityTools.jl are a collection of statistics that in some manner quantify the "association" between input datasets. Precisely what is meant by "association" depends on the measure, and precisely what is meant  by "quantify" depends on the *estimator* of that measure. 

Many association statistics exists. We have divided these statistics into (currently) three separate APIs, which are described below (use the left-hand menu to navigate the page).

- [Information measures](@ref)
- [Cross-map measures](@ref)
- [Correlation measures](@ref)

## [Overview](@id overview_association_measures)

| Type                    | Measure                               | Pairwise | Conditional | Function version               |
| ----------------------- | ------------------------------------- | :------: | :---------: | ------------------------------ |
| Correlation             | [`PearsonCorrelation`](@ref)          |    ✓    |     ✖      | [`pearson_correlation`](@ref)  |
| Correlation             | [`DistanceCorrelation`](@ref)         |    ✓    |     ✓      | [`distance_correlation`](@ref) |
| Closeness               | [`SMeasure`](@ref)                    |    ✓    |     ✖      | [`s_measure`](@ref)            |
| Closeness               | [`HMeasure`](@ref)                    |    ✓    |     ✖      | [`h_measure`](@ref)            |
| Closeness               | [`MMeasure`](@ref)                    |    ✓    |     ✖      | [`m_measure`](@ref)            |
| Closeness (ranks)       | [`LMeasure`](@ref)                    |    ✓    |     ✖      | [`l_measure`](@ref)            |
| Closeness               | [`JointDistanceDistribution`](@ref)   |    ✓    |     ✖      | [`jdd`](@ref)                  |
| Cross-mapping           | [`PairwiseAsymmetricInference`](@ref) |    ✓    |     ✖      | [`crossmap`](@ref)             |
| Cross-mapping           | [`ConvergentCrossMapping`](@ref)      |    ✓    |     ✖      | [`crossmap`](@ref)             |
| Conditional recurrence  | [`MCR`](@ref)                         |    ✓    |     ✖      | [`mcr`](@ref)                  |
| Conditional recurrence  | [`RMCD`](@ref)                        |    ✓    |     ✓      | [`rmcd`](@ref)                 |
| Shared information      | [`MIShannon`](@ref)                   |    ✓    |     ✖      | [`mutualinfo`](@ref)           |
| Shared information      | [`MIRenyiJizba`](@ref)                |    ✓    |     ✖      | [`mutualinfo`](@ref)           |
| Shared information      | [`MIRenyiSarbu`](@ref)                |    ✓    |     ✖      | [`mutualinfo`](@ref)           |
| Shared information      | [`MITsallisFuruichi`](@ref)           |    ✓    |     ✖      | [`mutualinfo`](@ref)           |
| Shared information      | [`PartialCorrelation`](@ref)          |    ✖    |     ✓      | [`partial_correlation`](@ref)  |
| Shared information      | [`CMIShannon`](@ref)                  |    ✖    |     ✓      | [`condmutualinfo`](@ref)       |
| Shared information      | [`CMIRenyiSarbu`](@ref)               |    ✖    |     ✓      | [`condmutualinfo`](@ref)       |
| Shared information      | [`CMIRenyiJizba`](@ref)               |    ✖    |     ✓      | [`condmutualinfo`](@ref)       |
| Information transfer    | [`TEShannon`](@ref)                   |    ✓    |     ✓      | [`transferentropy`](@ref)      |
| Information transfer    | [`TERenyiJizba`](@ref)                |    ✓    |     ✓      | [`transferentropy`](@ref)      |
| Part mutual information | [`PMI`](@ref)                         |    ✖    |     ✓      | [`pmi`](@ref)                  |
| Information asymmetry   | [`PA`](@ref)                          |    ✓    |     ✓      | [`asymmetry`](@ref)            |


## [Information measures](@id information_measures)

The information API defines bivariate and multivariate information measures, which we here define as *any functional of a multidimensional probability mass functions (PMFs) or multidimensional probability density*. Estimating an information measure from 
data is done by calling the [`information`](@ref) function with a relevant 
[`MultivariateInformationMeasure`](@ref) and a [`MultivariateInformationMeasureEstimator`](@ref). 
We have taken great care to make [`information`](@ref) modular. For one given measure, there may be multiple ways of estimating it. Individual docstrings detail estimation
possibilities.

```@docs
CausalityTools.information(::MultivariateInformationMeasureEstimator)
MultivariateInformationMeasure
MultivariateInformationMeasureEstimator
```

### Generic estimators

We provide a set of generic estimators that can be used to calculate 
potentially several types of information measures.

```@docs
JointProbabilities
EntropyDecomposition
MIDecomposition
CMIDecomposition
```


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
mutualinfo
MIShannon
MITsallisFuruichi
MITsallisMartin
MIRenyiJizba
MIRenyiSarbu
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

### Conditional mutual informations

```@docs
ConditionalMutualInformation
condmutualinfo
CMIShannon
CMIRenyiSarbu
CMIRenyiJizba
CMIRenyiPoczos
CMITsallis
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

### Partial mutual information

```@docs
PartialMutualInformation
```

### Transfer entropy

```@docs
TransferEntropy
transferentropy
TEShannon
TERenyiJizba
```


#### Transfer entropy estimators

```@docs
TransferEntropyEstimator
Zhu1
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

#### Counting and probabilities

For counting and probabilities, CausalityTools.jl extends the single-variable machinery
in ComplexityMeasures.jl to multiple variables.

```@docs
CausalityTools.Counts
CausalityTools.counts
CausalityTools.Probabilities
CausalityTools.probabilities
marginal
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

## [Correlation measures](@id correlation_measures)

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

## [Cross-map measures](@id cross_map_measures)

The cross-map measures define different ways of quantifying association based on the 
concept of "cross mapping", which has appeared in many contexts in the literature,
and gained huge popularity with  [Sugihara2012](@citet)'s on *convergent cross mapping*.

Since their paper, several cross mapping methods and frameworks have emerged in the
literature. In CausalityTools.jl, we provide a unified interface for using these cross
mapping methods.

To estimate a cross map measure, you simply input a [`CrossmapMeasure`](@ref) instance
as the first argument to a [`CrossmapEstimator`](@ref), which is then fed to 
the [`crossmap`](@ref) or [`predict`](@ref) functions. 

The Cross-map measures consists of the following functions.

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
