# Design philosophy

## Building blocks: probabilities and entropies

The building blocks of most information theoretic measures are [probabilities](@ref probabilities_header) and [entropies](@ref entropy_header). However, these quantities
also appear in many other applications, for example in complexity quantification.

Realizing this, we decided the gather the fundamentals required for estimation of
probabilities and entropies in a single package:
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl).
ComplexityMeasures.jl was built with modularity in mind, and provides a plethora of
estimators of probabilities and generalized entropies, both discrete and continuous.
These estimators are used frequently throughout CausalityTools.jl.

## The same names for different things: keeping track using *definitions*

In contrast to generalized entropies, which each have *one* definition, it gets a bit more
complicated when it comes to the "higher-level" measures we provide here.

Upon doing a literature review on the possible variants of information theoretic measures,
it become painstakingly obvious that authors use *the same name for different concepts*.
For beginners in the field of information theory and its practical applications, this
can be very confusing. In this package, we classify classify the various information
measures and their (sometimes non-equivalent) variants according to *definitions*.

This allows us to use Julia's multiple dispatch system to have a *single* function for
computing related, but disparate concepts, like [`mutualinfo`](@ref). As a bonus,
integration with probabilities estimators and entropy estimators from
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl)
is super simple, and makes it easy to develop and maintain the codebase.

### An example: mutual information

Quantities such as "mutual information" have multiple definitions in the
literature. For example, the Shannon mutual information (MI) has both a discrete and
continuous version, and for both cases, there are multiple equivalent mathematical formulas
for it. As an example, in this package, the [`MIShannon`](@ref) measure currently accepts 
two definitions. The
[`MIDefinitionShannonH3`](@ref) definition computes either discrete or continuous
MI by a sum of three entropy terms. Another example is the
[`MIDefinitionShannonDoubleSum`](@ref)), which computes discrete Shannon MI from
a double sum (for which joint probabilities are estimated using a contingency table).
Even more equivalent definitions are possible, for example by directly computing some form
of Kullback-Leibler divergence between joint distributions and products of marginal
distributions. Each such method comes with its own advantages/disadvantages when it
comes to practical estimation.

Moreover, Shannon MI is not the only mutual information! Many authors have proposed different
measures of shared information between random variables, taking as a starting point some
generalized entropy definition like [`Tsallis`](@ref). This has led to measures like
the "Tsallis mutual information". Like for Shannon MI, Tsallis MI also
has many definitions in the scientific literature, and *not all of them are equivalent*!
Naturally, this naming ambiguity has the potential to cause confusion when it comes to
the understanding and applicability of the methods, and the to interpretion of their results.

To alleviate potential confusion, and to streamline the estimation of these measures,
we have designed this package around *measures* and *definitions*. These measures 
are fed to a *single* function that unifies all related concepts. For mutual information,
this function is [`mutualinfo`](@ref).

For example,
we explicitly define a subtype of [`MutualInformationDefinition`](@ref)s for each
definition/variant of MI that appears in the literature. Specifying such a definition (e.g.
[`MIDefinitionShannonDoubleSum`](@ref)) to a particular measure (e.g. [`MIShannon`](@ref))
forces the computation does [`mutualinfo`](@ref) to follow that particular definition.

### What about other measures?

Mutual information is not the only measure that appear in multiple forms in the literature.
The same is the case for conditional mutual information (CMI), and for more basic
measures like the relative entropy. For practical purposes, this doesn't add any extra
complexity here. Each method (e.g. "conditional mutual information") has its
own function (e.g. [`condmutualinfo`](@ref)), which dispaches on some type of
conditional mutual information measure (e.g. [`CMIShannon`](@ref)), which can be
tweaked to use any [`ConditionalMutualInformationDefinition`].

For a complete overview of which estimators are compatible with which measures,
and to understand how computations are done, explore the docstrings for the definitions
and browse relevant sections in the menu!

### Universal estimation

Our modular package design comes with a super-convenient advantage: in principle,
one could estimate *any* information measure using *any* estimator. Although
the current interface doesn't allow *every* combination of measure
and estimator, you can already do a lot!

To compute your favorite measure, simply find a suitable estimator in one of the
overview tables, and apply it to some input data! Follow one of the plentiful
examples for inspiration.

If you're interested in a deeper understanding, we try to give mathematical formulas and implementation details as best we can in the docstrings of the various measures and
definitions.

### Summary

This package has been and is under heavy development. Don't hesitate to submit an
issue if you find something that doesn't work or doesn't make sense, or if there's
some functionality that you're missing!

Pull requests are also very welcome!
