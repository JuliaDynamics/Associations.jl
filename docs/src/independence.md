
## Independence testing

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
([`LocalPermutationTest`](@ref) test on your data with *over 20 different estimators* for
the conditional mutual information. If your application rather calls for the use
of traditional [surrogate data](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl),
the [`SurrogateTest`](@ref) test seamlessly integrates with any
time series surrogate method from TimeseriesSurrogates.jl.

### Independence test API

The independence test API is defined by

* [`independence`](@ref)
* [`IndependenceTest`](@ref)
* [`ConditionalIndependenceTest`](@ref)

### Independence tests

```@docs
independence
```

```@docs
ConditionalIndependenceTest
```

#### Surrogate test (global permutation)

```@docs
SurrogateTest
```

#### Local permutation

```@docs
LocalPermutationTest
```

### Examples

#### [`LocalPermutationTest`](@ref)

Here, we'll create a three-variable scenario where `X` and `Z` are connected through `Y`,
so that ``I(X; Z | Y) = 0`` and ``I(X; Y | Z) > 0``. We'll test for conditional
independence using Shannon conditional mutual information
([`CMIShannon`](@ref)). To estimate CMI, we'll use the [`Kraskov`](@ref) differential
entropy estimator, which naively computes CMI as a sum of entropy terms without guaranteed
bias cancellation.

```@example LOCAL_PERMUTATION_TEST
using CausalityTools

X = randn(1000)
Y = X .+ randn(1000) .* 0.4
Z = randn(1000) .+ Y
x, y, z = Dataset.((X, Y, Z))
test = LocalPermutationTest(CMIShannon(base = 2), Kraskov(k = 10), nsurr = 30)
test_result = independence(test, x, y, z)
```

We expect there to be a detectable influence from ``X`` to
``Y``, if we condition on ``Z`` or not, because ``Z`` doesn't influence neither ``X`` nor ``Y``.
The null hypothesis is that the first two variables are conditionally independent given the third, which we reject with a very low p-value. Hence, we accept the alternative
hypothesis that the first two variables ``X`` and ``Y``. are conditionally *dependent* given ``Z``.

```@example LOCAL_PERMUTATION_TEST
test_result = independence(test, x, z, y)
```

As expected, we cannot reject the null hypothesis that ``X`` and ``Z`` are conditionally independent given ``Y``, because ``Y`` is the variable that transmits information from
``X`` to ``Z``.

#### [`SurrogateTest`](@ref)

To demonstrate the [`SurrogateTest`](@ref) test, we use the transfer entropy measure,
which accepts either two input timeseries, or three timeseries when computing the
partial/conditional transfer entropy.

```@example surrogatecit_te
using CausalityTools
sys = logistic2_unidir(c_xy = 0.5) # x affects y, but not the other way around.
x, y = columns(trajectory(sys, 1000, Ttr = 10000))

test = SurrogateTest(TEShannon(), KSG1(k = 4))
independence(test, x, y)
```

As expected, we can reject the null hypothesis that the future of `y` is independent of 
`x`, because `x` does actually influence `y`. This doesn't change if we compute 
partial transfer entropy with respect to some random extra time series, because it
doesn't influence either variables.

```@example surrogatecit_te
independence(test, x, y, rand(length(x)))
```
