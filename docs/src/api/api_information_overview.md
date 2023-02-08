# [Information API](@id information_api)

This page outlines the information API. It contains a lot of information, so for
convenience, we list all concrete implementation of pairwise and conditional
association measures [here](@ref information_measures).

## [Design](@id information_measures_design)

We have taken great care to make sure that information estimators are reusable and modular.
Functions have the following general form.

```julia
f([measure], estimator, input_data...)

# Some examples
mutualinfo(MIShannon(base = ℯ), Kraskov(k = 1), x, y)
mutualinfo(MITsallisFuruichi(base = ℯ), KozachenkoLeonenko(k = 3), x, y)
condmutualinfo(CMIShannon(base = 2), ValueHistogram(3), x, y, z)
condmutualinfo(CMIRenyiJizba(base = 2), KSG2(k = 5), x, y, z)
condmutualinfo(CMIRenyiPoczos(base = 2), PoczosSchneiderCMI(k = 10), x, y, z)
```

This modular design really shines when it comes to independence testing and causal graph
inference. You can essentially test the performance of *any* independence `measure` with
*any* `estimator`, as long as their combination is implemented (and if it's not,
please submit a PR or issue!). We hope that this will both ease reproduction of
existing literature results, and spawn new research. Please let us know if you use the
package for something useful, or publish something based on it!

Information measures are either estimated using one of the following basic estimator types,

- [`ProbabilitiesEstimator`](@ref)s,
- [`DifferentialEntropyEstimator`](@ref)s,

or using measure-specific estimators:

- [`MutualInformationEstimator`](@ref)s are used with [`mutualinfo`](@ref)
- [`ConditionalMutualInformationEstimator`](@ref)s are used with [`condmutualinfo`](@ref)
- [`TransferEntropyEstimator`](@ref)s are used with [`transferentropy`](@ref)

## Naming convention: The same name for different things

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
