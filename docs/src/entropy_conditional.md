# Conditional entropy

## Conditional entropy API

The mutual information API is defined by

* [`ConditionalEntropy`](@ref),
* [`entropy_conditional`](@ref),

We provide a suite of estimators of various mutual information quantities. Many more
variants exist in the literature. Pull requests are welcome!

## Conditional entropy definitions

```@docs
ConditionalEntropy
CEShannon
CETsallisFuruichi
CETsallisAbe
```

More variants exist in the literature. Pull requests are welcome!

## Discrete conditional entropy

```@docs
entropy_conditional(::ConditionalEntropy, ::ContingencyMatrix)
```

### [Contingency matrix](@id contingency_matrix_ce)

Discrete conditional entropy can be computed directly from its sum-definition
by using the probabilities from a [`ContingencyMatrix`](@ref). This estimation
method works for  both numerical and categorical data, and the following
[`ConditionalEntropy`](@ref) definitions are supported.

|                             | [`ContingencyMatrix`](@ref) |
| --------------------------- | :-------------------------: |
| [`CEShannon`](@ref)         |             ✓              |
| [`CETsallisFuruichi`](@ref) |             ✓              |
| [`CETsallisAbe`](@ref)      |             ✓              |

### [Table of discrete conditional entropy estimators](@id probabilities_estimators_ce)

Here, we list the [`ProbabilitiesEstimator`](@ref)s that are compatible with
[`entropy_conditional`](@ref), and which definitions they are valid for.

| Estimator                    | Principle           | [`CEShannon`](@ref) | [`CETsallisAbe`](@ref) | [`CETsallisFuruichi`](@ref) |
| ---------------------------- | ------------------- | :-----------------: | :--------------------: | :-------------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         |         ✓          |           ✓           |              x              |
| [`ValueHistogram`](@ref)     | Binning (histogram) |         ✓          |           ✓           |              x              |
| [`SymbolicPermuation`](@ref) | Ordinal patterns    |         ✓          |           ✓           |              x              |
| [`Dispersion`](@ref)         | Dispersion patterns |         ✓          |           ✓           |              x              |

## Differential/continuous conditional entropy

### [Table of differential conditional entropy estimators](@id diffentropy_estimators_ce)

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

## Examples

### Discrete: example from Cover & Thomas

This is essentially example 2.2.1 in Cover & Thomas (2006), where they use the following
contingency table as an example. We'll take their example and manually construct
a [`ContingencyMatrix`](@ref) that we can use to compute the conditional entropy.
The [`ContingencyMatrix`](@ref) constructor takes the probabilities as the
first argument and the raw frequencies as the second argument.
Note also that Julia is column-major, so we need to transpose their example. Then their
`X` is in the first dimension of our contingency matrix (along columns) and their `Y` is
our second dimension (rows).

```@example ce_contingency_table
using CausalityTools
freqs_yx = [1//8 1//16 1//32 1//32; 
    1//16 1//8  1//32 1//32;
    1//16 1//16 1//16 1//16; 
    1//4  0//1  0//1  0//1];
freqs_xy = transpose(freqs_yx);
probs_xy = freqs_xy ./ sum(freqs_xy)
c_xy = ContingencyMatrix(probs_xy, freqs_xy)
```

The marginal distribution for `x` (first dimension) is

```@example ce_contingency_table
probabilities(c_xy, 1)
```

The marginal distribution for `y` (second dimension) is

```@example ce_contingency_table
probabilities(c_xy, 2)
```

And the Shannon conditional entropy ``H^S(X | Y)``

```@example ce_contingency_table
ce_x_given_y = entropy_conditional(CEShannon(), c_xy) |> Rational
```

This is the same as in their example. Hooray! To compute ``H^S(Y | X)``, we just need to
flip the contingency matrix.

```@example ce_contingency_table
probs_yx = freqs_yx ./ sum(freqs_yx);
c_yx = ContingencyMatrix(probs_yx, freqs_yx);
ce_y_given_x = entropy_conditional(CEShannon(), c_yx) |> Rational
```
