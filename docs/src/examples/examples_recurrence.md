# [Inferring directional influence using conditional recurrence](@id examples_recurrence)

## Computing the [`Recurrence`](@ref) measure for independent data

The interpretation of the [`Recurrence`](@ref) measure is that if two variables are
symmetrically coupled, then the conditional recurrence in both directions is equal.
Two variables that are uncoupled are symmetrically coupled (i.e. no coupling). We
therefore expect the difference in conditional recurrence to be around zero.

```@example
using CausalityTools
using StableRNGs
rng = StableRNG(1234)
x = rand(rng, 300)
y = rand(rng, 300)
m = Recurrence(r = 0.5)
Δ = conditional_recurrence(m, x, y) - conditional_recurrence(m, y, x)
```

This value is close to zero. To test if it is significantly indistinguishable from
zero, we can use a [`SurrogateAssociationTest`](@ref) (see example below).

## Independence test

```@example
using CausalityTools
using StableRNGs
rng = StableRNG(1234)
x = rand(rng, 300)
y = rand(rng, 300)
test = SurrogateAssociationTest(Recurrence(r = 0.5); rng, nshuffles = 100, surrogate = RandomShuffle())
independence(test, x, y)
```

As expected, we can't reject independence. What happens if two variables are coupled?

```@example
using CausalityTools
using StableRNGs
rng = StableRNG(1234)
x = rand(rng, 300)
z = x .+ rand(rng, 300)
test = SurrogateAssociationTest(Recurrence(r = 0.5); rng, nshuffles = 100, surrogate = RandomShuffle())
independence(test, x, z)
```

Now, because the variables are coupled, the evidence in the data support dependence.
