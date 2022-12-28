# Independence testing

```@docs
independence
```

```@docs
ConditionalIndependenceTest
```

## `LocalPermutation`

```@docs
LocalPermutation
```

### Examples

Here, we'll create a three-variable scenario where `X` and `Z` are connected through `Y`,
so that ``I(X; Z | Y) = 0`` and ``I(X; Y | Z) > 0``. We'll use the four-entropy-sum
definition ([`CMI4H`](@ref)) of Shannon conditional mutual information
([`CMIShannon`](@ref)), and use the [`Kraskov`](@ref) estimator to estimate each entropy
term.

```@example LOCAL_PERMUTATION_TEST
using CausalityTools
using Random

X = randn(10000)
Y = X .+ randn(10000)
Z = randn(10000) .+ 0.8*Y
x, y, z = Dataset.([X, Y, Z])

test = LocalPermutation(
    definition = CMI4H(),
    measure = CMIShannon(base = 2),
    est = Kraskov(k = 10))
test_result = independence(test, x, y, z)
```
