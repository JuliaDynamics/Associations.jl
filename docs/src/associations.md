# Association measures

The most basic components of CausalityTools.jl are a collection of statistics that in some manner quantify the "association" between input datasets. Precisely what is meant by "association" depends on the measure, and precisely what is meant by "quantify" depends on the *estimator* of that measure. 

To compute any association between two or more variable, use the [`association`](@ref) function with 
an [`AssociationMeasure`](@ref) (complete list in the table below) or an [`AssociationMeasureEstimator`](@ref).

## [Overview](@id overview_association_measures)

| Type                    | [`AssociationMeasure`](@ref)          | Pairwise | Conditional | Convenience function           |
|-------------------------|---------------------------------------|:--------:|:-----------:|--------------------------------|
| Correlation             | [`PearsonCorrelation`](@ref)          |    ✓     |      ✖      | [`pearson_correlation`](@ref)  |
| Correlation             | [`DistanceCorrelation`](@ref)         |    ✓     |      ✓      | [`distance_correlation`](@ref) |
| Closeness               | [`SMeasure`](@ref)                    |    ✓     |      ✖      | [`s_measure`](@ref)            |
| Closeness               | [`HMeasure`](@ref)                    |    ✓     |      ✖      | [`h_measure`](@ref)            |
| Closeness               | [`MMeasure`](@ref)                    |    ✓     |      ✖      | [`m_measure`](@ref)            |
| Closeness (ranks)       | [`LMeasure`](@ref)                    |    ✓     |      ✖      | [`l_measure`](@ref)            |
| Closeness               | [`JointDistanceDistribution`](@ref)   |    ✓     |      ✖      | [`jdd`](@ref)                  |
| Cross-mapping           | [`PairwiseAsymmetricInference`](@ref) |    ✓     |      ✖      | [`crossmap`](@ref)             |
| Cross-mapping           | [`ConvergentCrossMapping`](@ref)      |    ✓     |      ✖      | [`crossmap`](@ref)             |
| Conditional recurrence  | [`MCR`](@ref)                         |    ✓     |      ✖      | [`mcr`](@ref)                  |
| Conditional recurrence  | [`RMCD`](@ref)                        |    ✓     |      ✓      | [`rmcd`](@ref)                 |
| Shared information      | [`MIShannon`](@ref)                   |    ✓     |      ✖      | [`mutualinfo`](@ref)           |
| Shared information      | [`MIRenyiJizba`](@ref)                |    ✓     |      ✖      | [`mutualinfo`](@ref)           |
| Shared information      | [`MIRenyiSarbu`](@ref)                |    ✓     |      ✖      | [`mutualinfo`](@ref)           |
| Shared information      | [`MITsallisFuruichi`](@ref)           |    ✓     |      ✖      | [`mutualinfo`](@ref)           |
| Shared information      | [`PartialCorrelation`](@ref)          |    ✖     |      ✓      | [`partial_correlation`](@ref)  |
| Shared information      | [`CMIShannon`](@ref)                  |    ✖     |      ✓      | [`condmutualinfo`](@ref)       |
| Shared information      | [`CMIRenyiSarbu`](@ref)               |    ✖     |      ✓      | [`condmutualinfo`](@ref)       |
| Shared information      | [`CMIRenyiJizba`](@ref)               |    ✖     |      ✓      | [`condmutualinfo`](@ref)       |
| Information transfer    | [`TEShannon`](@ref)                   |    ✓     |      ✓      | [`transferentropy`](@ref)      |
| Information transfer    | [`TERenyiJizba`](@ref)                |    ✓     |      ✓      | [`transferentropy`](@ref)      |
| Part mutual information | [`PMI`](@ref)                         |    ✖     |      ✓      | [`pmi`](@ref)                  |
| Information asymmetry   | [`PA`](@ref)                          |    ✓     |      ✓      | [`asymmetry`](@ref)            |


```@docs
association
AssociationMeasure
AssociationMeasureEstimator
```

## [Information measures](@id information_api)

```@docs
MultivariateInformationMeasure
MultivariateInformationMeasureEstimator
```

### Generic information estimators

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
TEShannon
TERenyiJizba
```


#### Transfer entropy estimators

```@docs
TransferEntropyEstimator
Zhu1
Lindner
```


##### Convenience

```@docs
SymbolicTransferEntropy
Hilbert
Phase
Amplitude
```

##### Utilities

```@docs
optimize_marginals_te
EmbeddingTE
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

## [Correlation measures](@id correlation_api)

This page lists all available [`CorrelationMeasure`](@ref)s, as 
well as their convenience functions. The [examples](@ref correlation_examples)
is also useful.

### Pearson correlation

```@docs
PearsonCorrelation
```

### Partial correlation

```@docs
PartialCorrelation
```

### Distance correlation

```@docs
DistanceCorrelation
```

## [Cross-map measures](@id cross_map_api)

The cross-map measures define different ways of quantifying association based on the 
concept of "cross mapping", which has appeared in many contexts in the literature,
and gained huge popularity with  [Sugihara2012](@citet)'s on *convergent cross mapping*.

Since their paper, several cross mapping methods and frameworks have emerged in the
literature. In CausalityTools.jl, we provide a unified interface for using these cross
mapping methods.

To estimate a cross map measure, you simply input a [`CrossmapMeasure`](@ref) instance
as the first argument to a [`CrossmapEstimator`](@ref), which is then fed to 
the [`crossmap`](@ref) or [`predict`](@ref) functions. 

The cross-map measures consists of the following functions.

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

## [Closeness measures](@id closeness_api)

### Joint distance distribution

```@docs
JointDistanceDistribution
```

### S-measure

```@docs
SMeasure
```

### H-measure

```@docs
HMeasure
```

### M-measure

```@docs
MMeasure
```

### L-measure

```@docs
LMeasure
```

## Convenience

We encourage using [`association`](@ref) to compute associations between variables. However, since there 
are many types of established *names* for different types of associations, we provide the following 
convenience functions.

```@docs
joint_entropy
conditional_entropy
mutualinfo
condmutualinfo
transferentropy
pmi
pearson_correlation
partial_correlation
distance_correlation
mcr
rmcd
jdd
s_measure
h_measure
m_measure
l_measure
```