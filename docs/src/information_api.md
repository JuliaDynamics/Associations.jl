# [Information API](@id information_api)

This page outlines the information API. It contains a lot of information, so for
convenience, we list all concrete implementation of pairwise and conditional
association measures [here](@ref information_measures).

## [Design](@id information_measures_design)

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

## Estimators

Information measures are either estimated using one of the following basic estimator types,

- [`ProbabilitiesEstimator`](@ref)s,
- [`DifferentialEntropyEstimator`](@ref)s,

or using measure-specific estimators:

- [`MutualInformationEstimator`](@ref)s are used with [`mutualinfo`](@ref)
- [`ConditionalMutualInformationEstimator`](@ref)s are used with [`condmutualinfo`](@ref)
- [`TransferEntropyEstimator`](@ref)s are used with [`transferentropy`](@ref)

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
    distributions. Since these definitions are all equivalent, we only need once type
    ([`MIShannon`](@ref)) to represent them.
- But Shannon MI is not the  only type of mutual information! For example, "Tsallis mutual information"
    has been proposed in different variants by various authors. Despite sharing the
    same name, these are actually *nonequivalent definitions*. We've thus assigned
    them entirely different measure names (e.g. [`MITsallisFuruichi`](@ref) and
    [`MITsallisMartin`](@ref)), with the author name at the end.

## Probability mass functions (pmf)

Discrete information theoretic association measures and other quantities are estimated using
[`ProbabilitiesEstimator`](@ref)s. Here, we list probabilities estimators that are
compatible with this package.

| Estimator                                   | Principle                                      | Input data          |
| :------------------------------------------ | :--------------------------------------------- | :------------------ |
| [`Contingency`](@ref)                       | Count frequencies, optionally discretize first | `Any`               |
| [`CountOccurrences`](@ref)                  | Count of unique elements                       | `Any`               |
| [`ValueHistogram`](@ref)                    | Binning (histogram)                            | `Vector`, `Dataset` |
| [`TransferOperator`](@ref)                  | Binning (transfer operator)                    | `Vector`, `Dataset` |
| [`NaiveKernel`](@ref)                       | Kernel density estimation                      | `Dataset`           |
| [`SymbolicPermutation`](@ref)               | Ordinal patterns                               | `Vector`, `Dataset` |
| [`Dispersion`](@ref)                        | Dispersion patterns                            | `Vector`            |

### Contingency

```@docs
Contingency
```

### Count occurrences

```@docs
CountOccurrences
```

### Histograms (binning)

```@docs
ValueHistogram
RectangularBinning
FixedRectangularBinning
```

### Transfer operator (binning)

```@docs
TransferOperator
```

#### Utility methods/types

For explicit estimation of the transfer operator, see
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl).

```@docs
InvariantMeasure
invariantmeasure
transfermatrix
```

### Symbolic permutations

```@docs
SymbolicPermutation
```

### Dispersion patterns

```@docs
Dispersion
```

### Kernel density

```@docs
NaiveKernel
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

### Probabilities API

The probabilities API is defined by

- [`ProbabilitiesEstimator`](@ref)
- [`probabilities`](@ref)
- [`probabilities_and_outcomes`](@ref)
- [`ContingencyMatrix`](@ref)
- [`contingency_matrix`](@ref)

#### Probabilities and outcomes

```@docs
ProbabilitiesEstimator
probabilities
probabilities!
Probabilities
probabilities_and_outcomes
outcomes
outcome_space
total_outcomes
missing_outcomes
```

#### Encodings

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

##### Available encodings

```@docs
OrdinalPatternEncoding
GaussianCDFEncoding
RectangularBinEncoding
```

## Contingency tables

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

### Contingency matrix API

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

## Conditional entropy API

The conditional entropy API is defined by

- [`ConditionalEntropy`](@ref),
- [`entropy_conditional`](@ref),

```@docs
ConditionalEntropy
entropy_conditional
```

### Conditional entropy definitions

```@docs
CEShannon
CETsallisFuruichi
CETsallisAbe
```

## [Mutual information API](@id api_mutualinfo)

The mutual information API is defined by

- [`MutualInformation`](@ref),
- [`mutualinfo`](@ref),
- [`MutualInformationEstimator`](@ref).

We provide a suite of estimators of various mutual information quantities. Many more
variants exist in the literature. Pull requests are welcome!

```@docs
MutualInformation
mutualinfo
```

### Mutual information estimators

```@docs
MutualInformationEstimator
KraskovStögbauerGrassberger1
KraskovStögbauerGrassberger2
GaoKannanOhViswanath
GaoOhViswanath
```

## Conditional mutual information API

The condition mutual information (CMI) API is defined by

- [`ConditionalMutualInformation`](@ref),
- [`mutualinfo`](@ref),
- [`ConditionalMutualInformationEstimator`](@ref).

### CMI definitions

```@docs
ConditionalMutualInformation
condmutualinfo
```

### Dedicated CMI estimators

```@docs
ConditionalMutualInformationEstimator
FPVP
MesnerShalisi
PoczosSchneiderCMI
Rahimzamani
```

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
