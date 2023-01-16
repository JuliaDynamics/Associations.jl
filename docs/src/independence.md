# Independence testing

A common application of information theoretic methods such as conditional mutual
information ([`condmutualinfo`](@ref)) is in the context of null hypothesis testing
for the conditional independence of variables.

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
of traditional [surrogate data](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl),
the [`SurrogateCIT`](@ref) test seamlessly integrates with any
time series surrogate method from TimeseriesSurrogates.jl.

## Independence test API

The independence test API is defined by

* [`conditional_independence`](@ref)
* [`ConditionalIndependenceTest`](@ref)

## Independence tests

```@docs
conditional_independence
```

```@docs
ConditionalIndependenceTest
```

### Surrogate test (global permutation)

```@docs
SurrogateCIT
```

### Local permutation

```@docs
LocalPermutation
```

## Examples

### [`LocalPermutation`](@ref)

Here, we'll create a three-variable scenario where `X` and `Z` are connected through `Y`,
so that ``I(X; Z | Y) = 0`` and ``I(X; Y | Z) > 0``. We'll test for conditional
independence using Shannon conditional mutual information
([`CMIShannon`](@ref)). To estimate CMI, we'll use the [`Kraskov`](@ref) differential
entropy estimator, which naively computes CMI as a sum of entropy terms without guaranteed
bias cancellation.

```@example LOCAL_PERMUTATION_TEST
using CausalityTools

X = randn(2000)
Y = X .+ randn(2000) .* 0.4
Z = randn(2000) .+ Y
x, y, z = Dataset.((X, Y, Z))
test = LocalPermutation(nsurr = 35,
    measure = CMIShannon(base = 2),
    est = Kraskov(k = 10)
);
test_result = conditional_independence(test, x, y, z)
```

We expect there to be a detectable influence from ``X`` to
``Y`` even conditioning on ``Z``, because ``Z`` doesn't influence neither ``X`` nor ``Y``.
The null hypothesis is that the first two variables are conditionally independent given the third, which we reject with a very low p-value. Hence, we accept the alternative
hypothesis that the first two variables ``X`` and ``Y``. are conditionally *dependent* given ``Z``.

```@example LOCAL_PERMUTATION_TEST
test_result = conditional_independence(test, x, z, y)
```

As expected, we cannot reject the null hypothesis that ``X`` and ``Z`` are conditionally independent given ``Y``, because ``Y`` is the variable that transmits information from
``X`` to ``Z``.
