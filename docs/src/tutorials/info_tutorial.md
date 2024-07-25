# [Information measure tutorial](@id tutorial_infomeasures)

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
the *definition* of the measure, which is then used as input to an 
*estimator* of that measure. This estimator is then given to [`association`](@ref).

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

# Specify a discretization. We discretize per column.
disc = CodifyVariables(OrdinalPatterns(m=2))

def = KLDivergence()
est = JointProbabilities(def, disc)
association(est, x, y) # should be close to 0
```

Divergences are examples of *asymmetric* information measures, which we can see by 
flipping the order of the input data.

```@example INFO_TUTORIAL
association(est, y, x)
```

## Conditional entropies

[Conditional entropies](@ref conditional_entropies) are another example of asymmetric
information measures. They all have in common that 
they are functions of a joint pmf, and can therefore also be estimated using the
[`JointProbabilities`](@ref) estimator. This time, we'll use a rectangular binning
with 3 bins along each dimension to discretize the data.

```@example INFO_TUTORIAL
using CausalityTools
using Random; rng = Xoshiro(1234)

x, y = randn(rng, 1000), randn(rng, 1000)
disc = CodifyVariables(ValueBinning(3))
def = ConditionalEntropyShannon(base = 2)
est = JointProbabilities(def, disc)
association(est, x, y)
```

## Joint entropies

[Joint entropies](@ref joint_entropies), on the other hand, are *symmetric*. Joint
entropies are functionals of a joint pmf, so we can still use the
[`JointProbabilities`](@ref) estimator. This time, we use a [`Dispersion`](@ref)
based discretization.

```@example INFO_TUTORIAL
using CausalityTools
using Random; rng = Xoshiro(1234)

x, y = randn(rng, 1000), randn(rng, 1000)
disc = CodifyVariables(Dispersion())
est = JointProbabilities(JointEntropyShannon(base = 2), disc)
association(est, x, y) ≈ association(est, y, x) # should be true
```

## Mutual informations

Mutual informations, in particular [`MIShannon`](@ref) is an often-used symmetric 
measure for quantifing the (possibly nonlinear) association between variables. It appears
in both  discrete and differential form, and can be estimated in a multitude of ways. For 
example, one can use dedicated [`MutualInformationEstimator`](@ref)s such as 
[`KraskovStögbauerGrassberger2`](@ref) or [`GaussianMI`](@ref):

```@example INFO_TUTORIAL
using CausalityTools
using StateSpaceSets
# We'll construct two state space sets, to illustrate how we can discretize 
# multidimensional inputs using `CodifyPoints`.
x, y = StateSpaceSet(rand(rng, 1000, 2)), StateSpaceSet(rand(rng, 1000, 2))
est = KSG1(MIShannon(base = 2), k = 10)
association(est, x, y)
```

The result should be symmetric:

```@example INFO_TUTORIAL
association(est, x, y) ≈ association(est, y, x) # should be true
```

One can also estimate mutual information using the [`EntropyDecomposition`](@ref) 
estimator, or (like above) using the [`JointProbabilities`](@ref) estimator.
Let's construct a differential entropy based estimator based on the [`Kraskov`](@ref)
estimator.

```@example INFO_TUTORIAL
est_diff = EntropyDecomposition(MIShannon(base = 2), Kraskov(Shannon(), k=10))
```

```@example INFO_TUTORIAL
association(est_diff, x, y)
```

We can also construct a discrete entropy based estimator based on e.g. [`PlugIn`](@ref)
estimator of [`Shannon`](@ref) entropy.

```@example INFO_TUTORIAL
# We know that `x` and `y` were generated from a uniform distribution above,
# so we set the minimum and maximum values of the encoding to 0 and 1,
# respectively.
encoding = RelativeMeanEncoding(0.0, 1.0; n = 4)
disc = CodifyPoints(encoding)
est_disc = JointProbabilities(MIShannon(base = 2), disc)
association(est_disc, x, y)
```
