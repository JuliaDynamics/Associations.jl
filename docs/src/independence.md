# Independence testing

A common application of information theoretic methods such as conditional mutual
information ([`condmutualinfo`](@ref)) is in the context of null hypothesis testing
for the conditional independence of variables. In time series applications,
[`transferentropy`](@ref), which is just a variant of ([`condmutualinfo`](@ref)),
the same applies!

Depending on the context, the input data and the method used, there are many considerations
to be made about how to perform this conditional independence testing. Luckily, many
excellent frameworks for doing so exist in the literature.

Here, we present some commonly used independence tests from the scientific literature,
which can all be seamlessly used with the function [`independence`](@ref),
with *any* measure that quantifies conditional independence, in combination
with *any* compatible estimator. 

For example, in just a few lines of code, you can perform Runge's local permutation 
([`LocalPermutation`](@ref) test on your data with *over 20 different estimators* for
the conditional mutual information. If your application rather calls for the use
of [surrogate data](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl),
the [`GlobalPermutation`](@ref) test seamlessly integrates with any
time series surrogate method from TimeseriesSurrogates.jl.

## API

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
