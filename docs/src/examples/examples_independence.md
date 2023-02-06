# [Independence testing](@id examples_independence)

## [`LocalPermutationTest`](@ref)

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

## [`SurrogateTest`](@ref)

### [[`TEShannon`](@ref)](@id examples_surrogatetest_teshannon)

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

### [[`SMeasure`](@ref)](@id examples_surrogatetest_smeasure)

```@example quickstart_smeasure
using CausalityTools
x, y = randn(3000), randn(3000)
measure = SMeasure(dx = 3, dy = 3)
s = s_measure(measure, x, y)
```

The `s` statistic is larger when there is stronger coupling and smaller
when there is weaker coupling. To check whether `s` is significant (i.e. large
enough to claim directional dependence), we can use a [`SurrogateTest`](@ref),
like [here](@ref examples_surrogatetest_smeasure).

```@example quickstart_smeasure
test = SurrogateTest(measure)
independence(test, x, y)
```

The p-value is high, and we can't reject the null at any reasonable significance level.
Hence, there isn't evidence in the data to support directional coupling from `x` to `y`.

What happens if we use coupled variables?

```@example quickstart_smeasure
z = x .+ 0.1y
independence(test, x, z)
```

Now we can confidently reject the null (independence), and conclude that there is
evidence in the data to support directional dependence from `x` to `z`.
