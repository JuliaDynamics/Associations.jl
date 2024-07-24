# Association measures

## Associations API

The most basic components of CausalityTools.jl are a collection of statistics that in some manner quantify the "association" between input datasets. Precisely what is meant by "association" depends on the measure, and precisely what is meant by "quantify" depends on the *estimator* of that measure. We formalize this notion below with the [`association`](@ref)
function, which dispatches on [`AssociationMeasureEstimator`](@ref) and [`AssociationMeasure`](@ref).


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
CMITsallisPapapetrou
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

### Encoding API

A fundamental operation when computing multivariate information measures from data is *discretization*.  The following functions and types are used by CausalityTools.jl to perform discretization of input data.

```@docs
codify
Discretization
```

#### Encoding per variable/column

```@docs
CodifyVariables
```

The sliding-window discretization is formally done by applying some [`OutcomeSpace`](@ref) to each variable/column. Pick between the following outcome spaces

```@docs
OutcomeSpace
UniqueElements
CosineSimilarityBinning
Dispersion
OrdinalPatterns
BubbleSortSwaps
ValueBinning
RectangularBinning
FixedRectangularBinning
```

#### Encoding per sample/row

```@docs
CodifyPoints
```

```@docs
Encoding
GaussianCDFEncoding
OrdinalPatternEncoding
RelativeMeanEncoding
RelativeFirstDifferenceEncoding
UniqueElementsEncoding
RectangularBinEncoding
CombinationEncoding
```

### [Counts and probabilities](@id counts_and_probabilities_api)

For counting and probabilities, CausalityTools.jl extends the single-variable machinery
in ComplexityMeasures.jl to multiple variables.

```@docs
CausalityTools.Counts
CausalityTools.counts(::OutcomeSpace)
```

```@docs
CausalityTools.Probabilities
CausalityTools.probabilities(::OutcomeSpace)
```

#### Utility functions

```@docs
marginal
```

### Single-variable information API

Below we list some relevant types from
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl) that 
are used for the [`EntropyDecomposition`](@ref) estimator.

#### Single-variable information measures

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
Ensemble
```

### Advanced utility methods

For most use cases, it is sufficient to provide a [`CrossmapEstimator`](@ref) to 
[`association`](@ref) to compute a cross map measure. However, in some cases it 
can be useful to have more fine-grained controls. We offer a few utility functions
for this purpose.

In the example where we [reproduce Figures 3C and 3D](@ref example_sugihara_figs3Cand3D) of [Sugihara2012](@ref), these lower-level 
functions are used.

```@docs
predict
crossmap
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

## Recurrence measures

```@docs
MCR
RMCD
```