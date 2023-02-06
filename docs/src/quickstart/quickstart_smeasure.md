# [S-measure](@id quickstart_smeasure)

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

```@example quickstart_jdd
test = SurrogateTest(measure)
independence(test, x, y)
```

The p-value is high, and we can't reject the null at any reasonable significance level.

```@example
```
