# [Information API](@id information_api)

This page outlines the information API. It contains a lot of information, so for
convenience, we list all concrete implementation of pairwise and conditional
association measures [here](@ref information_measures).

## [API & design](@id information_measures)

### Modularity

We have taken great care to make sure that information estimators are reusable and modular.
Functions have the following general form.

```julia
f([measure], estimator, input_data...)

# Some examples
mutualinfo(MIShannon(base = ℯ), Kraskov(k = 1), x, y)
mutualinfo(MITsallisFuruichi(base = ℯ), KozachenkoLeonenko(k = 3), x, y)
condmutualinfo(CMIShannon(base = 2), ValueHistogram(3), x, y, z)
condmutualinfo(CMIRenyiJizba(base = 2), KSG2(k = 5), x, y, z)
condmutualinfo(CMIRenyiPoczos(base = 2), PoczosSchneiderCMI(k = 10), x, y, z)
```

This modular design really shines when it comes to independence testing and causal graph
inference. You can essentially test the performance of *any* independence `measure` with
*any* `estimator`, as long as their combination is implemented (and if it's not,
please submit a PR or issue!). We hope that this will both ease reproduction of
existing literature results, and spawn new research. Please let us know if you use the
package for something useful, or publish something based on it!

### Naming convention: The same name for different things

Upon doing a literature review on the possible variants of information theoretic measures,
it become painstakingly obvious that authors use *the same name for different concepts*.
For novices, and experienced practitioners too, this can be confusing.
Our API clearly distinguishes between methods that are conceptually the same but named
differently in the literature due to differing *estimation* strategies, from methods
that actually have different definitions.

- Multiple, equivalent definitions occur for example for the Shannon mutual
    information (MI; [`MIShannon`](@ref)), which has both a discrete and continuous version, and there there are multiple equivalent mathematical formulas for them: a direct sum/integral
    over a joint probability mass function (pmf), as a sum of three entropy terms, and as
    a Kullback-Leibler divergence between the joint pmf and the product of the marginal
    distributions. Since these are all equivalent, we only need once type (`[`MIShannon`](@ref)) to represent them.
- But Shannon MI is not the  only type of mutual information! "Tsallis mutual information"
    has been proposed in different variants by various authors. Despite sharing the
    same name, these are actually *nonequivalent definitions*. Naming ambiguities like
    these are likely to cause confusion. We've thus assigned
    them entirely different measure names (e.g. [`MITsallisFuruichi`](@ref) and
    [`MITsallisMartin`](@ref)), with the author name at the end.

Every measure starts with an abbrevation of the quantity it measures, followed by the name of the measure,
for example:

- [`MIShannon`](@ref) has many *equivalent* definitions that all share the same name.
- [`MITsallisFuruichi`](@ref) and [`MITsallisMartin`](@ref) are separate measures,
    because they are defined by *nonequivalent* mathematical formulas.

## Probability mass functions (pmf)

### Probabilities API

The probabilities API is defined by

- [`ProbabilitiesEstimator`](@ref)
- [`probabilities`](@ref)
- [`probabilities_and_outcomes`](@ref)
- [`ContingencyMatrix`](@ref)
- [`contingency_matrix`](@ref)

and related functions that you will find in the following documentation blocks:

#### Probabilitities

```@docs
ProbabilitiesEstimator
probabilities
probabilities!
Probabilities
```

#### Outcomes

```@docs
probabilities_and_outcomes
outcomes
outcome_space
total_outcomes
missing_outcomes
```

### Estimators

#### [Overview of probabilities estimators](@id probabilities_estimators)

Any of the following estimators can be used with [`probabilities`](@ref)
(in the column "input data"  it is assumed that the `eltype` of the input is `<: Real`).
Some estimators can also be used with [`contingency_matrix`](@ref) to estimate
multivariate pmfs.

| Estimator                                   | Principle                                      | Input data          |
| :------------------------------------------ | :--------------------------------------------- | :------------------ |
| [`Contingency`](@ref)                       | Count frequencies, optionally discretize first | `Any`               |
| [`CountOccurrences`](@ref)                  | Count of unique elements                       | `Any`               |
| [`ValueHistogram`](@ref)                    | Binning (histogram)                            | `Vector`, `Dataset` |
| [`TransferOperator`](@ref)                  | Binning (transfer operator)                    | `Vector`, `Dataset` |
| [`NaiveKernel`](@ref)                       | Kernel density estimation                      | `Dataset`           |
| [`SymbolicPermutation`](@ref)               | Ordinal patterns                               | `Vector`, `Dataset` |
| [`SymbolicWeightedPermutation`](@ref)       | Ordinal patterns                               | `Vector`, `Dataset` |
| [`SymbolicAmplitudeAwarePermutation`](@ref) | Ordinal patterns                               | `Vector`, `Dataset` |
| [`SpatialSymbolicPermutation`](@ref)        | Ordinal patterns in space                      | `Array`             |
| [`Dispersion`](@ref)                        | Dispersion patterns                            | `Vector`            |
| [`SpatialDispersion`](@ref)                 | Dispersion patterns in space                   | `Array`             |
| [`Diversity`](@ref)                         | Cosine similarity                              | `Vector`            |
| [`WaveletOverlap`](@ref)                    | Wavelet transform                              | `Vector`            |
| [`PowerSpectrum`](@ref)                     | Fourier transform                              | `Vector`            |

### Contingency

```@docs
Contingency
```

### Count occurrences

```@docs
CountOccurrences
```

### Histograms

```@docs
ValueHistogram
RectangularBinning
FixedRectangularBinning
```

### Symbolic permutations

```@docs
SymbolicPermutation
SymbolicWeightedPermutation
SymbolicAmplitudeAwarePermutation
```

### Dispersion patterns

```@docs
Dispersion
```

### Transfer operator (binning)

```@docs
TransferOperator
```

For explicit estimation of the transfer operator, see
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl).

#### Utility methods/types

```@docs
InvariantMeasure
invariantmeasure
transfermatrix
```

### Kernel density

```@docs
NaiveKernel
```

### Local likelihood

```@docs
LocalLikelihood
```

### Timescales

```@docs
WaveletOverlap
PowerSpectrum
```

### Diversity

```@docs
Diversity
```

### Spatial estimators

```@docs
SpatialSymbolicPermutation
SpatialDispersion
```

### Encodings

#### Encodings API

Some probability estimators first "encode" input data into an intermediate representation indexed by the positive integers. This intermediate representation is called an "encoding".

The encodings API is defined by:

- [`Encoding`](@ref)
- [`encode`](@ref)
- [`decode`](@ref)

```@docs
Encoding
encode
decode
```

#### Available encodings

```@docs
OrdinalPatternEncoding
GaussianCDFEncoding
RectangularBinEncoding
```

### Contingency tables

To estimate discrete information theoretic quantities that are functions of more than
one variable, we must estimate empirical joint probability mass functions (pmf).
The function [`contingency_matrix`](@ref) accepts an arbitrary number of equal-length
input data and returns the corresponding multidimensional contingency table as a
[`ContingencyMatrix`](@ref). From this table, we can extract the necessary joint and
marginal pmfs for computing any discrete function of multivariate discrete probability
distributions. This is essentially the multivariate analogue of
[`Probabilities`](@ref).

But why would I use a [`ContingencyMatrix`](@ref) instead of some other indirect estimation
method, you may ask. The answer is that [`ContingencyMatrix`](@ref) allows you to
compute *any* of the information theoretic quantities offered in this package for *any*
type of input data. You input data can literally be any hashable type, for example `String`,
`Tuple{Int, String, Int}`, or `YourCustomHashableDataType`.

In the case of numeric data, using a [`ContingencyMatrix`](@ref) is typically a
bit slower than other dedicated estimation procedures.
For example, quantities like discrete Shannon-type [`condmutualinfo`](@ref) are faster to
estimate using a formulation based on sums of four entropies (the H4-principle). This
is faster because we can both utilize the blazingly fast [`Dataset`](@ref) structure directly,
and we can avoid *explicitly* estimating the entire joint pmf, which demands many
extra calculation steps. Whatever you use in practice depends on your use case and
available estimation methods, but you can always fall back to contingency matrices
for any discrete measure.

#### Contingency matrix API

```@docs
ContingencyMatrix
contingency_matrix
marginal_encodings
```

## [Entropies](@id entropies)

### Entropies API

The entropies API is defined by

- [`EntropyDefinition`](@ref)
- [`entropy`](@ref)
- [`DifferentialEntropyEstimator`](@ref)

Please be sure you have read the [Terminology](@ref) section before going through the API here, to have a good idea of the different "flavors" of entropies and how they all come together over the common interface of the [`entropy`](@ref) function.

### Entropy definitions

```@docs
EntropyDefinition
Shannon
Renyi
Tsallis
Kaniadakis
Curado
StretchedExponential
```

### Discrete

```@docs
entropy(::EntropyDefinition, ::ProbabilitiesEstimator, ::Any)
entropy_maximum
entropy_normalized
```

### Differential/continuous

```@docs
entropy(::EntropyDefinition, ::DifferentialEntropyEstimator, ::Any)
```

#### Table of differential entropy estimators

The following estimators are *differential* entropy estimators, and can also be used
with [`entropy`](@ref).

Each [`DifferentialEntropyEstimator`](@ref)s uses a specialized technique to approximate relevant
densities/integrals, and is often tailored to one or a few types of generalized entropy.
For example, [`Kraskov`](@ref) estimates the [`Shannon`](@ref) entropy.

| Estimator                    | Principle         | Input data | [`Shannon`](@ref) | [`Renyi`](@ref) | [`Tsallis`](@ref) | [`Kaniadakis`](@ref) | [`Curado`](@ref) | [`StretchedExponential`](@ref) |
| :--------------------------- | :---------------- | :--------- | :---------------: | :-------------: | :---------------: | :------------------: | :--------------: | :----------------------------: |
| [`KozachenkoLeonenko`](@ref) | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Kraskov`](@ref)            | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Zhu`](@ref)                | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`ZhuSingh`](@ref)           | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Gao`](@ref)                | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Goria`](@ref)              | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Lord`](@ref)               | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Vasicek`](@ref)            | Order statistics  | `Vector`   |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Ebrahimi`](@ref)           | Order statistics  | `Vector`   |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Correa`](@ref)             | Order statistics  | `Vector`   |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`AlizadehArghami`](@ref)    | Order statistics  | `Vector`   |        ✓         |        x        |         x         |          x           |        x         |               x                |

```@docs
DifferentialEntropyEstimator
```

```@docs
Kraskov
KozachenkoLeonenko
Zhu
ZhuSingh
Gao
Goria
Lord
Vasicek
AlizadehArghami
Ebrahimi
Correa
```


## Conditional entropy

### Conditional entropy API

The conditional entropy API is defined by

- [`ConditionalEntropy`](@ref),
- [`entropy_conditional`](@ref),

### Conditional entropy definitions

```@docs
ConditionalEntropy
CEShannon
CETsallisFuruichi
CETsallisAbe
```

More variants exist in the literature. Pull requests are welcome!

### Discrete conditional entropy

```@docs
entropy_conditional(::ConditionalEntropy, ::ContingencyMatrix)
```

#### [Contingency matrix](@id contingency_matrix_ce)

Discrete conditional entropy can be computed directly from its sum-definition
by using the probabilities from a [`ContingencyMatrix`](@ref). This estimation
method works for  both numerical and categorical data, and the following
[`ConditionalEntropy`](@ref) definitions are supported.

|                             | [`ContingencyMatrix`](@ref) |
| --------------------------- | :-------------------------: |
| [`CEShannon`](@ref)         |             ✓              |
| [`CETsallisFuruichi`](@ref) |             ✓              |
| [`CETsallisAbe`](@ref)      |             ✓              |

#### [Table of discrete conditional entropy estimators](@id probabilities_estimators_ce)

Here, we list the [`ProbabilitiesEstimator`](@ref)s that are compatible with
[`entropy_conditional`](@ref), and which definitions they are valid for.

| Estimator                    | Principle           | [`CEShannon`](@ref) | [`CETsallisAbe`](@ref) | [`CETsallisFuruichi`](@ref) |
| ---------------------------- | ------------------- | :-----------------: | :--------------------: | :-------------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         |         ✓          |           ✓           |              x              |
| [`ValueHistogram`](@ref)     | Binning (histogram) |         ✓          |           ✓           |              x              |
| [`SymbolicPermuation`](@ref) | Ordinal patterns    |         ✓          |           ✓           |              x              |
| [`Dispersion`](@ref)         | Dispersion patterns |         ✓          |           ✓           |              x              |

### Differential/continuous conditional entropy

#### [Table of differential conditional entropy estimators](@id diffentropy_estimators_ce)

Continuous/differential mutual information may be estimated using any of our
[`DifferentialEntropyEstimator`](@ref)s that support multivariate input data.

| Estimator                        | Principle         | [`CEShannon`](@ref) | [`CETsallisAbe`](@ref) | [`CETsallisFuruichi`](@ref) |
| -------------------------------- | ----------------- | :-----------------: | :--------------------: | :-------------------------: |
| [`Kraskov`](@ref)                | Nearest neighbors |         ✓          |           x           |              x              |
| [`Zhu`](@ref)                    | Nearest neighbors |         ✓          |           x           |              x              |
| [`ZhuSingh`](@ref)               | Nearest neighbors |         ✓          |           x           |              x              |
| [`Gao`](@ref)                    | Nearest neighbors |         ✓          |           x           |              x              |
| [`Goria`](@ref)                  | Nearest neighbors |         ✓          |           x           |              x              |
| [`Lord`](@ref)                   | Nearest neighbors |         ✓          |           x           |              x              |
| [`LeonenkoProzantoSavani`](@ref) | Nearest neighbors |         ✓          |           x           |              x              |

## Mutual information

### [Mutual information API](@id api_mutualinfo)

The mutual information API is defined by

- [`MutualInformation`](@ref),
- [`mutualinfo`](@ref),
- [`MutualInformationEstimator`](@ref).

We provide a suite of estimators of various mutual information quantities. Many more
variants exist in the literature. Pull requests are welcome!

#### Definitions

```@docs
MutualInformation
```

#### Dedicated estimators

```@docs
mutualinfo(est::MutualInformationEstimator, ::Any, ::Any) 
```

```@docs
MutualInformationEstimator
KraskovStögbauerGrassberger1
KraskovStögbauerGrassberger2
GaoKannanOhViswanath
GaoOhViswanath
```

##### [Table of dedicated estimators](@id dedicated_estimators_mi)

| Estimator                              |    Type    |     Principle     | [`MIShannon`](@ref) | [`MITsallisFuruichi`](@ref) | [`MITsallisMartin`](@ref) | [`MIRenyiSarbu`](@ref) | [`MIRenyiJizba`](@ref) |
| -------------------------------------- | :--------: | :---------------: | :-----------------: | :-------------------------: | :-----------------------: | :--------------------: | :--------------------: |
| [`KraskovStögbauerGrassberger1`](@ref) | Continuous | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`KraskovStögbauerGrassberger2`](@ref) | Continuous | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`GaoKannanOhViswanath`](@ref)         |   Mixed    | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`GaoOhViswanath`](@ref)               | Continuous | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |

#### Discrete mutual information

```@docs
mutualinfo(::ProbabilitiesEstimator, ::Any, ::Any)
```

##### [Table of discrete mutual information estimators](@id @id dedicated_probabilities_estimators_mi)

Here, we list the [`ProbabilitiesEstimator`](@ref)s that can be used to compute discrete
[`mutualinformation`](@ref).

| Estimator                    | Principle           | [`MIShannon`](@ref) | [`MITsallisFuruichi`](@ref) | [`MITsallisMartin`](@ref) | [`MIRenyiJizba`](@ref) | [`MIRenyiSarbu`](@ref) |
| ---------------------------- | ------------------- | :-----------------: | :-------------------------: | :-----------------------: | :--------------------: | :--------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         |         ✓          |             ✓              |            ✓             |           ✓           |           x           |
| [`ValueHistogram`](@ref)     | Binning (histogram) |         ✓          |             ✓              |            ✓             |           ✓           |           x           |
| [`SymbolicPermuation`](@ref) | Ordinal patterns    |         ✓          |             ✓              |            ✓             |           ✓           |           x           |
| [`Dispersion`](@ref)         | Dispersion patterns |         ✓          |             ✓              |            ✓             |           ✓           |           x           |

##### [Contingency matrix](@id contingency_matrix_mi)

```@docs
mutualinfo(::MutualInformation, ::ContingencyMatrix)
```

Discrete mutual information can be computed directly from its double-sum definition
by using the probabilities from a [`ContingencyMatrix`](@ref). This estimation
method works for    both numerical and categorical data, and the following
[`MutualInformation`](@ref)s are supported.

|                             | [`ContingencyMatrix`](@ref) |
| --------------------------- | :-------------------------: |
| [`MIShannon`](@ref)         |             ✓              |
| [`MITsallisFuruichi`](@ref) |             ✓              |
| [`MITsallisMartin`](@ref)   |             ✓              |
| [`MIRenyiSarbu`](@ref)      |             ✓              |
| [`MIRenyiJizba`](@ref)      |             ✓              |

#### Differential/continuous mutual information

```@docs
mutualinfo(::DifferentialEntropyEstimator, ::Any, ::Any)
```

##### [Table of differential mutual information estimators](@id dedicated_diffentropy_estimators_mi)

In addition to the dedicated differential mutual information estimators listed above,
continuous/differential mutual information may also be estimated using any of our
[`DifferentialEntropyEstimator`](@ref) that support multivariate input data.
When using these estimators, mutual information is computed as a sum
of entropy terms (with different dimensions), and no bias correction is applied.

| Estimator                        | Principle         | [`MIShannon`](@ref) | [`MITsallisFuruichi`](@ref) | [`MITsallisMartin`](@ref) | [`MIRenyiJizba`](@ref) | [`MIRenyiSurbu`](@ref) |
| -------------------------------- | ----------------- | :-----------------: | :-------------------------: | :-----------------------: | :--------------------: | :--------------------: |
| [`Kraskov`](@ref)                | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Zhu`](@ref)                    | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`ZhuSingh`](@ref)               | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Gao`](@ref)                    | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Goria`](@ref)                  | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Lord`](@ref)                   | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`LeonenkoProzantoSavani`](@ref) | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |

## Conditional mutual information

### CMI API

The condition mutual information API is defined by

- [`ConditionalMutualInformation`](@ref),
- [`mutualinfo`](@ref),
- [`ConditionalMutualInformationEstimator`](@ref).

### CMI definitions

```@docs
ConditionalMutualInformation
```

### Dedicated CMI estimators

```@docs
condmutualinfo(::ConditionalMutualInformationEstimator, ::Any, ::Any, ::Any)
```

```@docs
ConditionalMutualInformationEstimator
FPVP
MesnerShalisi
PoczosSchneiderCMI
Rahimzamani
```

#### [Table of dedicated CMI estimators](@id condmutualinfo_dedicated_estimators)

| Estimator                    | Principle         | [`CMIShannon`](@ref) | [`CMIRenyiPoczos`](@ref) |
| ---------------------------- | ----------------- | :------------------: | :----------------------: |
| [`FPVP`](@ref)               | Nearest neighbors |          ✓          |            x             |
| [`MesnerShalisi`](@ref)      | Nearest neighbors |          ✓          |            x             |
| [`Rahimzamani`](@ref)        | Nearest neighbors |          ✓          |            x             |
| [`PoczosSchneiderCMI`](@ref) | Nearest neighbors |          x           |            ✓            |
| [`GaussianCMI`](@ref)        | Parametric        |          ✓          |            x             |

### Estimation through mutual information

```@docs
condmutualinfo(::MutualInformationEstimator, ::Any, ::Any, ::Any)
```

| Estimator                              |    Type    |     Principle     | [`CMIShannon`](@ref) |
| -------------------------------------- | :--------: | :---------------: | :------------------: |
| [`KraskovStögbauerGrassberger1`](@ref) | Continuous | Nearest neighbors |          ✓          |
| [`KraskovStögbauerGrassberger2`](@ref) | Continuous | Nearest neighbors |          ✓          |
| [`GaoKannanOhViswanath`](@ref)         |   Mixed    | Nearest neighbors |          ✓          |
| [`GaoOhViswanath`](@ref)               | Continuous | Nearest neighbors |          ✓          |
| [`GaussianMI`](@ref)                   |            |    Parametric     |          ✓          |

### Discrete CMI

```@docs
condmutualinfo(::ProbabilitiesEstimator, ::Any, ::Any, ::Any)
```

#### [Table of discrete mutual information estimators](@id mutualinfo_overview)

Here, we list the [`ProbabilitiesEstimator`](@ref)s that are compatible with
[`condmutualinfo`](@ref), and which definitions they are valid for.

| Estimator                    | Principle           | [`CMIShannon`](@ref) | [`CMIRenyiSarbu`](@ref) |
| ---------------------------- | ------------------- | :------------------: | :---------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         |          ✓          |           ✓            |
| [`ValueHistogram`](@ref)     | Binning (histogram) |          ✓          |           ✓            |
| [`SymbolicPermuation`](@ref) | Ordinal patterns    |          ✓          |           ✓            |
| [`Dispersion`](@ref)         | Dispersion patterns |          ✓          |           ✓            |

### Differential CMI

```@docs
condmutualinfo(::DifferentialEntropyEstimator, ::Any, ::Any, ::Any)
```

| Estimator                        | Principle         | Input data | [`CMIShannon`](@ref) |
| -------------------------------- | ----------------- | ---------- | :------------------: |
| [`Kraskov`](@ref)                | Nearest neighbors | `Dataset`  |          ✓          |
| [`Zhu`](@ref)                    | Nearest neighbors | `Dataset`  |          ✓          |
| [`Gao`](@ref)                    | Nearest neighbors | `Dataset`  |          ✓          |
| [`Goria`](@ref)                  | Nearest neighbors | `Dataset`  |          ✓          |
| [`Lord`](@ref)                   | Nearest neighbors | `Dataset`  |          ✓          |
| [`LeonenkoProzantoSavani`](@ref) | Nearest neighbors | `Dataset`  |          ✓          |

## Transfer entropy

The transfer entropy API is made up of the following functions and types,
which are listed below:

- [`transferentropy`](@ref).
- [`TransferEntropy`](@ref), and its subtypes.
- [`EmbeddingTE`](@ref), which exists to provide embedding instructions to
    subtypes of [`TransferEntropy`](@ref).
- [`TransferEntropyEstimator`](@ref), and its subtypes.

### Transfer entropy API

```@docs
transferentropy
EmbeddingTE
optimize_marginals_te
```

#### Transfer entropy definitions

```@docs
TransferEntropy
```

#### Transfer entropy estimators

```@docs
TransferEntropyEstimator
Zhu1
Lindner
```

#### Convenience

```@docs
SymbolicTransferEntropy
```

```@docs
Hilbert
Phase
Amplitude
```
