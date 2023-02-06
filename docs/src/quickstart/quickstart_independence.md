# [Independence testing](@id quickstart_independence)

## [[`JointDistanceDistributionTest`](@ref)](@id quickstart_jddtest)

Let's use the built-in `logistic2_bidir` discrete dynamical system to create a pair of
bidirectionally coupled time series and use the [`JointDistanceDistributionTest`](@ref)
to see if we can confirm from observed time series that these variables are
bidirectionally coupled. We'll use a significance level of `1 - α = 0.99`, i.e. `α = 0.01`.

We start by generating some time series and configuring the test.

```@example quickstart_jddtest_logistic
using CausalityTools
sys = logistic2_bidir(c_xy = 0.5, c_yx = 0.4)
x, y = columns(trajectory(sys, 2000, Ttr = 10000))
measure = JointDistanceDistribution(D = 5, B = 5)
test = JointDistanceDistributionTest(measure)
```

Now, we test for independence in both directions.

```@example quickstart_jddtest_logistic
independence(test, x, y)
```

```@example quickstart_jddtest_logistic
independence(test, y, x)
```

As expected, the null hypothesis is rejected in both directions at the pre-determined 
significance level, and hence we detect directional coupling in both directions. But what happens if there is no coupling?

```@example quickstart_jddtest_logistic
sys = logistic2_bidir(c_xy = 0.00, c_yx = 0.0)
x, y = columns(trajectory(sys, 1000, Ttr = 10000));
rxy = independence(test, x, y)
ryx = independence(test, y, x)
pvalue(rxy), pvalue(ryx)
```

At significance level `0.99`, we can't reject the null in either direction, hence there's not
enough evidence in the data to suggest directional coupling.

## [[`SurrogateTest`](@ref)](@id quickstart_surrogatetest)

### Mutual information

#### Categorical

In this example, we expect the `preference` and the `food` variables to be independent.

```@example
using CausalityTools
# Simulate 
n = 1000
preference = rand(["yes", "no"], n)
food = rand(["veggies", "meat", "fish"], n)
test = SurrogateTest(MIShannon(), Contingency())
independence(test, preference, food)
```

As expected, there's not enough evidence to reject the null hypothesis that the
variables are independent.

### Conditional mutual information

#### Categorical

Here, we simulate a survey at a ski resort. The data are such that the place a person
grew up is associated with how many times they fell while going skiing. The control
happens through an intermediate variable `preferred_equipment`, which indicates what
type of physical activity the person has engaged with in the past. Some activities
like skateboarding leads to better overall balance, so people that are good on
a skateboard also don't fall, and people that to less challenging activities fall
more often.

We should be able to reject `places ⫫ experience`, but not reject
`places ⫫ experience | preferred_equipment`.  Let's see if we can detect these
relationships using (conditional) mutual information.

```@example indep_cmi
using CausalityTools
n = 10000

places = rand(["city", "countryside", "under a rock"], n);
preferred_equipment = map(places) do place
    if cmp(place, "city") == 1
        return rand(["skateboard", "bmx bike"])
    elseif cmp(place, "countryside") == 1
        return rand(["sled", "snowcarpet"])
    else
        return rand(["private jet", "ferris wheel"])
    end
end;
experience = map(preferred_equipment) do equipment
    if equipment ∈ ["skateboard", "bmx bike"]
        return "didn't fall"
    elseif equipment ∈ ["sled", "snowcarpet"]
        return "fell 3 times or less"
    else
        return "fell uncontably many times"
    end
end;

test_mi = independence(SurrogateTest(MIShannon(), Contingency()), places, experience)
```

As expected, the evidence favors the alternative hypothesis that `places` and 
`experience` are dependent.

```@example  indep_cmi
test_cmi = independence(SurrogateTest(CMIShannon(), Contingency()), places, experience, preferred_equipment)
```

Again, as expected, when conditioning on the mediating variable, the dependence disappears,
and we can't reject the null hypothesis of independence.

## [[`LocalPermutationTest`](@ref)](@id quickstart_localpermutationtest)

### [`CMIShannon`](@ref)

#### Continuous data

```@example
using CausalityTools
x = randn(1000)
y = randn(1000) .+ 0.7 .* x
z = sin.(randn(1000)) .* 0.5 .* y
test = LocalPermutationTest(CMIShannon(; base = 2), FPVP(k = 10), nshuffles = 30)
independence(test, x, z, y)
```
