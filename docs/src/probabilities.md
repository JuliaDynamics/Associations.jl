# Probability mass functions (pmf)

## API

The probabilities API is defined by

- [`ProbabilitiesEstimator`](@ref)
- [`probabilities`](@ref)
- [`probabilities_and_outcomes`](@ref)
- [`ContingencyMatrix`](@ref)
- [`contingency_matrix`](@ref)

and related functions that you will find in the following documentation blocks:

### Probabilitities

```@docs
ProbabilitiesEstimator
probabilities
probabilities!
Probabilities
```

### Outcomes

```@docs
probabilities_and_outcomes
outcomes
outcome_space
total_outcomes
missing_outcomes
```

## Estimators

### [Overview of probabilities estimators](@id probabilities_estimators)

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

### TransferOperator (binning)

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

## Encodings

### Encodings API

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

### Available encodings

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

### Examples

The following two mutual information estimates on the integer ("categorical") vectors `x`, `y` are equivalent (up to rounding errors), but the latter is much faster to compute.

```@example
using CausalityTools
n = 1000
x = rand(1:3, n)
y = rand(1:4, n)
mcont = mutualinfo(MIShannon(), contingency_matrix(x, y))
mfast = mutualinfo(MIShannon(), CountOccurrences(), x, y)
mcont, mfast
```

On purely categorical data, you *have* to use the contingency matrix approach.

```julia
using CausalityTools
n = 100
likeit = rand(["yes", "no"], n)
food = rand(["veggies", "meat", "fish"], n)
mutualinfo(MIShannon(), contingency_matrix(likeit, food))
```
