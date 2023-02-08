# Closeness measures

## [S-measure](@id quickstart_smeasure)

### Computing the `s`-statistic

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

## [Joint distance distribution](@id quickstart_jdd)

### Computing the `Δ`-statistic

```@example quickstart_jdd
using CausalityTools
x, y = randn(3000), randn(3000)
measure = JointDistanceDistribution(D = 3, B = 5)
Δ = jdd(measure, x, y)
```

The joint distance distribution measure indicates directional coupling between
`x` and `y` if `Δ` is skewed towards positive values. We can use a [`JointDistanceDistributionTest`](@ref) to formally check this.

```@example quickstart_jdd
test = JointDistanceDistributionTest(measure)
independence(test, x, y)
```

The p-value is fairly low, and depending on the significance level `1 - α`, we cannot
reject the null hypothesis that `Δ` is not skewed towards positive values, and hence
we cannot reject that the variables are independent.