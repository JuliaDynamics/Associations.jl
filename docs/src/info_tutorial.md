# [Information measure tutorial](@id info_tutorial)

CausalityTools.jl extends the single-variate information API in
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl)
to information measures of multiple variables. 

## Definitions

We define **"information measure"** as some functional of probability 
mass functions or probability densities. This definition may or may not agree with
literature usage, depending on the context. We made this choice pragmatically based on
user-friendlyness and coding-friendlyness, but still trying to maintain some
level of meaningful terminology.

## Basic strategy

To *estimate* a multivariate information measure in practice, you must first specify
the [definition](@ref) of the measure, which is then used as input to an 
[estimator](@ref). This estimator is then given to [`information`](@ref), or one 
of the [convenience methods](@ref convenience_info).

!!! note "Naming convention: The same name for different things"
    Upon doing a literature review on the possible variants of information theoretic measures,
    it become painstakingly obvious that authors use *the same name for different concepts*.
    For novices, and experienced practitioners too, this can be confusing.
    Our API clearly distinguishes between methods that are conceptually the same but named
    differently in the literature due to differing *estimation* strategies, from methods
    that actually have different definitions.

    - Multiple, equivalent definitions occur for example for the Shannon mutual
        information (MI; [`MIShannon`](@ref)), which has both a discrete and continuous version, and there there are multiple equivalent mathematical formulas for them: a direct sum/integral
        over a joint probability mass function (pmf), as a sum of three entropy terms, and as
        a Kullback-Leibler divergence between the joint pmf and the product of the marginal
        distributions. Since these definitions are all equivalent, we only need once type
        ([`MIShannon`](@ref)) to represent them.
    - But Shannon MI is not the  only type of mutual information! For example, "Tsallis mutual information"
        has been proposed in different variants by various authors. Despite sharing the
        same name, these are actually *nonequivalent definitions*. We've thus assigned
        them entirely different measure names (e.g. [`MITsallisFuruichi`](@ref) and
        [`MITsallisMartin`](@ref)), with the author name at the end.


## Distances/divergences

There are many information measures in the literature that aim to quantify the 
distance/divergence between two probability mass functions (pmf) or densities. You can 
find those that we implement [here](@ref divergences_and_distances).

As an example, let's quantify the [`KLDivergence`](@ref) between two probability 
mass functions estimated by symbolizing two input vectors `x` and `y` using 
[`OrdinalPatterns`](@ref). Since the discrete [`KLDivergence`](@ref) can be 
expressed as a function of a joint pmf, we can use the [`JointProbabilities`](@ref)
estimator.

```@example INFO_TUTORIAL
using CausalityTools
using Random; rng = MersenneTwister(1234)
x, y = rand(rng, 1000), rand(rng, 1000)
est = JointProbabilities(KLDivergence(), OrdinalPatterns(m=2))
information(est, x, y) # should be close to 0
```

Divergences are examples of *asymmetric* information measures, which we can see by 
flipping the order of the input data.

```@example INFO_TUTORIAL
information(est, y, x)
```

## Conditional entropies

[Conditional entropies](@ref conditional_entropies) are another example of asymmetric
information measures. They all have in common that 
they are functions of a joint pmf, and can therefore also be estimated using the
[`JointProbabilities`](@ref) estimator. This time, we'll use a rectangular binning
with 3 bins along each dimension to discretize the data.

```@example INFO_TUTORIAL
x, y = randn(rng, 1000), randn(rng, 1000)
est = JointProbabilities(ConditionalEntropyShannon(base = 2), ValueBinning(3))
information(est, x, y)
```

## Joint entropies

[Joint entropies](@ref joint_entropies), on the other hand, are *symmetric*. Joint
entropies are functionals of a joint pmf, so we can still use the
[`JointProbabilities`](@ref) estimator. This time, we use a [`Dispersion`](@ref)
based discretization.

```@example INFO_TUTORIAL
x, y = randn(rng, 1000), randn(rng, 1000)
est = JointProbabilities(JointEntropyShannon(base = 2), Dispersion())
information(est, x, y) == information(est, y, x) # should be true
```

## Mutual informations

Mutual informations, in particular [`MIShannon`](@ref) is an often-used symmetric 
measure for quantifing the (possibly nonlinear) association between variables. It appears
in both  discrete and differential form, and can be estimated in a multitude of ways. For 
example, one can use dedicated [`MutualInformationEstimator`](@ref)s such as 
[`KSG2`](@ref) or [`GaussianMI`](@ref):

```@example INFO_TUTORIAL
x, y = randn(rng, 1000), randn(rng, 1000)
est = KSG1(MIShannon(base = 2), k = 10)
information(est, x, y)
```

The result should be symmetric:

```@example INFO_TUTORIAl
information(est, x, y) == information(est, y, x) # should be true
```

One can also estimate mutual information using the [`EntropyDecomposition`](@ref) 
estimator, or (like above) using the [`JointProbabilities`](@ref) estimator.
Let's construct a differential entropy based estimator based on the [`Kraskov`](@ref)
estimator.

```@example INFO_TUTORIAL
est_diff = EntropyDecomposition(MIShannon(base = 2), Kraskov(Shannon(), k=10))
```

We can also construct a discrete entropy based estimator based on e.g. [`PlugIn`](@ref)
estimator of [`Shannon`](@ref) entropy.

```@example INFO_TUTORIAL
est_disc = EntropyDecomposition(MIShannon(base = 2), PlugIn(Shannon()), ValueBinning(2))
```

These estimators use different estimation methods, so give different results:

```@example INFO_TUTORIAL
information(est_diff, x, y)
```

```@example INFO_TUTORIAL
information(est_disc, x, y)
```
