# Conditional entropy

The function that estimates conditional entropy from data is [`entropy_conditional`](@ref).
It can estimate different types of conditional entropies from the scientific literature,
and each is represented here as a subtype of [`ConditionalEntropy`](@ref). Because a
conditional entropy can be formulated in many different ways, each conditional entropy
type can be estimated according to multiple definitiosn, which are represented by
subtypes of [`ConditionalEntropyDefinition`](@ref)).

To see which estimators are compatible with the various definitions, see the
[overview table](@ref conditionalentropy_overview) below.

## Conditional entropy API

```@docs
entropy_conditional
ConditionalEntropy
ConditionalEntropyDefinition
```

## Definitions

More variants of the discrete Tsallis conditional entropy are possible beyond those
defined below.
Pull requests are welcome!

### Shannon conditional entropy

```@docs
CEShannon
```

### Tsallis conditional entropy

```@docs
CETsallisFuruichi
CETsallisAbe
```

## Discrete

### [Table of discrete conditional entropy estimators] (@id conditionalentropy_overview)

Here, we list the [`ProbabilitiesEstimator`](@ref)s that are compatible with
[`entropy_conditional`](@ref), and which definitions they are valid for.

| Estimator                    | Principle           |     Input data     | [`CEShannon`](@ref) | [`CETsallisFuruichi`](@ref) | [`CETsallisAbe`](@ref) |
| ---------------------------- | ------------------- | :----------------: | :-----------------: | :-------------------------: | :--------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         | `Vector`/`Dataset` |         ✅          |             ✅              |           ✅           |
| [`ValueHistogram`](@ref)     | Binning (histogram) | `Vector`/`Dataset` |         ✅          |             ✅              |           ✅           |
| [`SymbolicPermuation`](@ref) | Ordinal patterns    |      `Vector`      |         ✅          |             ✅              |           ✅           |
| [`Dispersion`](@ref)         | Dispersion patterns |      `Vector`      |         ✅          |             ✅              |           ✅           |

## Differential/continuous

We currently don't support the computation of differential conditional entropy.
Pull requests are welcome!

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
ce_y_given_x = entropy_conditional(CEShannon(), c_yx) |> Rational
```

This is the same as in their example. Hooray! To compute ``H^S(Y | X)``, we just need to
flip the contingency matrix.

```@example ce_contingency_table
probs_yx = freqs_yx ./ sum(freqs_yx);
c_yx = ContingencyMatrix(probs_yx, freqs_yx);
ce_x_given_x = entropy_conditional(CEShannon(), c_xy) |> Rational
```
