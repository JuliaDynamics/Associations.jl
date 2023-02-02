
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
([`LocalPermutationTest`](@ref) test on your data with *over 20 different estimators* for
the conditional mutual information. If your application rather calls for the use
of traditional [surrogate data](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl),
the [`SurrogateTest`](@ref) test seamlessly integrates with any
time series surrogate method from TimeseriesSurrogates.jl.

## Independence test API

The independence test API is defined by

* [`independence`](@ref)
* [`IndependenceTest`](@ref)
* [`ConditionalIndependenceTest`](@ref)

```@docs
independence
```

## Independence tests

```@docs
ConditionalIndependenceTest
```

### Surrogate test (global permutation)

```@docs
SurrogateTest
```

### Local permutation

```@docs
LocalPermutationTest
```
