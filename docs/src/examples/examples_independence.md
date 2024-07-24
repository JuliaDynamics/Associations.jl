# Examples of independence testing


## [[`CorrTest`](@ref)](@id example_CorrTest)

```@example corrtest_example
using CausalityTools
using Random; rng = StableRNG(1234)

# Some normally distributed data
X = randn(rng, 1000) 
Y = 0.5*randn(rng, 1000) .+ X
Z = 0.5*randn(rng, 1000) .+ Y
W = randn(rng, 1000);
```

Let's test a few independence relationships. For example, we expect that `X ⫫ W`.
We also expect dependence `X !⫫ Z`, but this dependence should vanish when
conditioning on the intermediate variable, so we expect `X ⫫ Z | Y`.

```@example corrtest_example
independence(CorrTest(), X, W)
```

As expected, the outcome is that we can't reject the null hypothesis that `X ⫫ W`.

```@example corrtest_example
independence(CorrTest(), X, Z)
```

However, we *can* reject the  null hypothesis that `X ⫫ Z`, so the evidence favors
the alternative hypothesis `X !⫫ Z`.

```@example corrtest_example
independence(CorrTest(), X, Z, Y)
```

As expected, the correlation between `X` and `Z` significantly vanishes when conditioning
on `Y`, because `Y` is solely responsible for the observed correlation between `X` and `Y`.

## [[`JointDistanceDistributionTest`](@ref)](@id examples_independence_JointDistanceDistributionTest)

### Bidirectionally coupled logistic maps

Let's use the built-in `logistic2_bidir` discrete dynamical system to create a pair of
bidirectionally coupled time series and use the [`JointDistanceDistributionTest`](@ref)
to see if we can confirm from observed time series that these variables are
bidirectionally coupled. We'll use a significance level of `1 - α = 0.99`, i.e. `α = 0.01`.

```@example examples_independence_JointDistanceDistributionTest
using CausalityTools
using DynamicalSystemsBase
using Random; rng = Xoshiro(1234)
Base.@kwdef struct Logistic2Bidir{V, C1, C2, R1, R2, Σx, Σy, R}
    xi::V = [0.5, 0.5]
    c_xy::C1 = 0.1
    c_yx::C2 = 0.1
    r₁::R1 = 3.78
    r₂::R2 = 3.66
    σ_xy::Σx = 0.05
    σ_yx::Σy = 0.05
    rng::R = Random.default_rng()
end

function system(definition::Logistic2Bidir)
    return DiscreteDynamicalSystem(eom_logistic2bidir, definition.xi, definition)
end

function eom_logistic2bidir(u, p::Logistic2Bidir, t)
    (; xi, c_xy, c_yx, r₁, r₂, σ_xy, σ_yx, rng) = p
    x, y = u
    f_xy = (y +  c_xy*(x + σ_xy * rand(rng)) ) / (1 + c_xy*(1+σ_xy))
    f_yx = (x +  c_yx*(y + σ_yx * rand(rng)) ) / (1 + c_yx*(1+σ_yx))
    dx = r₁ * (f_yx) * (1 - f_yx)
    dy = r₂ * (f_xy) * (1 - f_xy)
    return SVector{2}(dx, dy)
end
```

We start by generating some time series and configuring the test.

```@example examples_independence_JointDistanceDistributionTest
using CausalityTools
sys = system(Logistic2Bidir(c_xy = 0.5, c_yx = 0.4))
x, y = columns(first(trajectory(sys, 2000, Ttr = 10000)))
measure = JointDistanceDistribution(D = 5, B = 5)
test = JointDistanceDistributionTest(measure)
```

Now, we test for independence in both directions.

```@example examples_independence_JointDistanceDistributionTest
independence(test, x, y)
```

```@example examples_independence_JointDistanceDistributionTest
independence(test, y, x)
```

As expected, the null hypothesis is rejected in both directions at the pre-determined 
significance level, and hence we detect directional coupling in both directions.

### Non-coupled logistic maps

What happens in the example above if there is no coupling?

```@example examples_independence_JointDistanceDistributionTest
sys = system(Logistic2Bidir(c_xy = 0.00, c_yx = 0.0))
x, y = columns(first(trajectory(sys, 1000, Ttr = 10000)));
rxy = independence(test, x, y)
ryx = independence(test, y, x)
pvalue(rxy), pvalue(ryx)
```

At significance level `0.99`, we can't reject the null in either direction, hence there's not
enough evidence in the data to suggest directional coupling.

## [[`SurrogateAssociationTest`](@ref)](@id examples_surrogatetest)


### [Distance correlation](@id example_SurrogateAssociationTest_DistanceCorrelation)

```@example
using CausalityTools
x = randn(1000)
y = randn(1000) .+ 0.5x
independence(SurrogateAssociationTest(DistanceCorrelation()), x, y)
```

### [Partial correlation](@id example_SurrogateAssociationTest_PartialCorrelation)

```@example
using CausalityTools
x = randn(1000)
y = randn(1000) .+ 0.5x
z = randn(1000) .+ 0.8y
independence(SurrogateAssociationTest(PartialCorrelation()), x, z, y)
```


### [[`SMeasure`](@ref)](@id example_SurrogateAssociationTest_SMeasure)

```@example example_SurrogateAssociationTest_SMeasure
using CausalityTools
x, y = randn(1000), randn(1000)
measure = SMeasure(dx = 4, dy = 3)
s = association(measure,     x, y)
```

The `s` statistic is larger when there is stronger coupling and smaller
when there is weaker coupling. To check whether `s` is significant (i.e. large
enough to claim directional dependence), we can use a [`SurrogateAssociationTest`](@ref).

```@example example_SurrogateAssociationTest_SMeasure
test = SurrogateAssociationTest(measure)
independence(test, x, y)
```

The p-value is high, and we can't reject the null at any reasonable significance level.
Hence, there isn't evidence in the data to support directional coupling from `x` to `y`.

What happens if we use coupled variables?

```@example example_SurrogateAssociationTest_SMeasure
z = x .+ 0.1y
independence(test, x, z)
```

Now we can confidently reject the null (independence), and conclude that there is
evidence in the data to support directional dependence from `x` to `z`.


### [[`MIShannon`](@ref), categorical](@id example_SurrogateAssociationTest_MIShannon_categorical)

In this example, we expect the `preference` and the `food` variables to be independent.

```@example example_SurrogateAssociationTest_MIShannon_categorical
using CausalityTools
using Random; rng = Xoshiro(1234)
# Simulate 
n = 1000
preference = rand(rng, ["yes", "no"], n)
food = rand(rng, ["veggies", "meat", "fish"], n)
est = JointProbabilities(MIShannon(), CodifyVariables(UniqueElements()))
test = SurrogateAssociationTest(est)
independence(test, preference, food)
```

As expected, there's not enough evidence to reject the null hypothesis that the
variables are independent.

### [[`CMIShannon`](@ref), categorical](@id example_SurrogateAssociationTest_CMIShannon_categorical)

Here, we simulate a survey at a ski resort. The data are such that the place a person
grew up is associated with how many times they fell while going skiing. The control
happens through an intermediate variable `preferred_equipment`, which indicates what
type of physical activity the person has engaged with in the past. Some activities
like skateboarding leads to better overall balance, so people that are good on
a skateboard also don't fall, and people that to less challenging activities fall
more often.

We should be able to reject `places ⫫ experience`, but not reject
`places ⫫ experience | preferred_equipment`.  Let's see if we can detect these
relationships using (conditional) mutual information.

```@example indep_cmi
using CausalityTools
using Random; rng = Xoshiro(1234)
n = 1000

places = rand(rng, ["city", "countryside", "under a rock"], n);
preferred_equipment = map(places) do place
    if cmp(place, "city") == 1
        return rand(rng, ["skateboard", "bmx bike"])
    elseif cmp(place, "countryside") == 1
        return rand(rng, ["sled", "snowcarpet"])
    else
        return rand(rng, ["private jet", "ferris wheel"])
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
est_mi = JointProbabilities(MIShannon(), CodifyVariables(UniqueElements()))
test = SurrogateAssociationTest(est_mi)
independence(test, places, experience)
```

As expected, the evidence favors the alternative hypothesis that `places` and 
`experience` are dependent.

```@example  indep_cmi
est_cmi = JointProbabilities(CMIShannon(), CodifyVariables(UniqueElements()))
test = SurrogateAssociationTest(est_cmi)
independence(test, places, experience, preferred_equipment)
```

Again, as expected, when conditioning on the mediating variable, the dependence disappears,
and we can't reject the null hypothesis of independence.
