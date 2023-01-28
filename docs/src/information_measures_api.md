
# [API & design](@id information_measures)

Information measures are build on [probabilities](@ref probabilities_header)/densities
and [entropies](@ref entropy_header). We implement estimators of these quantities in
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl).
ComplexityMeasures.jl was built with modularity in mind, and provides a plethora of
estimators of probabilities and generalized entropies, both discrete and continuous.
These estimators are used frequently throughout CausalityTools.jl, relyin on the fact
that any "high-level" information measure, in some way or another, can be expressed
in terms of probabilities or entropies.

- Information measures are computed in their discrete form by using
    [`ProbabilitiesEstimator`](@ref).
- Information measures are computed in their differential/continuous
    form by using [`DifferentialEntropyEstimator`](@ref)s. Many measures also
    have dedicated estimators (like [`MutualInformationEstimator`](@ref) for
    [`mutualinfo`](@ref), some of which are designed to compute continuous quantities.

## Naming convention: The same name for different things

In contrast to generalized entropies, which each have *one* definition, it gets a bit more
complicated when it comes to the "higher-level" measures we provide here.

Upon doing a literature review on the possible variants of information theoretic measures,
it become painstakingly obvious that authors use *the same name for different concepts*.
For novices in the field of information theory, this can be very confusing. This package
is designed to alleviate any confusion regarding the names of information theoretic
quantities and their estimation.

We first consider the straight-forward case of multiple definitions. The Shannon mutual
information (MI) has both a discrete and continuous version, and
there there are multiple equivalent mathematical formulas for them: a direct sum/integral
over a joint probability mass function (pmf), as a sum of three entropy terms, and as
a Kullback-Leibler divergence between the joint pmf and the product of the marginal
distributions. Since these are all equivalent, we only need once type (`[`MIShannon`](@ref))
to represent them.

But Shannon MI variant this is not the only type of mutual information! Taking as a
starting point some generalized entropy definition like [`Tsallis`](@ref), several authors
have proposed variants of "Tsallis mutual information". Like Shannon MI, the Tsallis
variant also has many definitions in the scientific literature, and *not all of them are
equivalent, even though they are referred to by the same name*! Naming ambiguities like
these are likely to cause confusion.

To alleviate any confusion, we group *equivalent* definitions of a information measure
in a single type. Nonequivalent definitions are assigned separate types. Every measure
starts with an abbrevation of the quantity it measures, followed by the name of the
measure: [`CERenyi`](@ref) measures conditional Rényi entropy, and [`CEShannon`](@ref)
measures conditional Shannon entropy. If there are multiple definitions for the
same name, the author name is appended to the type:
[`MITsallisFuruichi`](@ref) and [`MITsallisMartin`](@ref) are separate measures,
because they are defined by *nonequivalent* mathematical formulas, whereas
[`MIShannon`](@ref) has many *equivalent* definitions.

To estimate some information measure, an instance of the measure type (e.g.
[`MIShannon`](@ref)) is combined with an *estimator*, which control *how* the quantity
is computed, given som input data. The most basic estimators are
[`ProbabilitiesEstimator`](@ref)s for discrete measures, and
[`DifferentialEntropyEstimator`](@ref)s for continuous/differential measures.
Some measures have dedicated estimators, that may be discrete, continuous or try to
estimate a mixture of discrete and continuous data.

### Examples

[Here](@ref example_mi_quickstart))'s an example of computing Shannon mutual information
using various estimators on various kinds of data.

Other measure like [`condmutualinfo`](@ref) also have multiple estimation routes.
To compute your favorite measure, simply find a suitable estimator in one of the
overview tables, and apply it to some input data! Follow one of the examples for
inspiration.

### Summary

With this modular API, one could in principle estimate *any* information measure using *any* estimator. Although the current interface doesn't allow *every* combination of measure
and estimator (and it's probably not theoretically meaningful to do so),
you can already do a lot!

If you're interested in a deeper understanding, we try to give mathematical formulas and
implementation details as best we can in the docstrings of the various measures and
definitions.
