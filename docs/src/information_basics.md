# Design philosophy

## Building blocks

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
can be very confusing. This package is designed to alleviate any confusion regarding the names of information theoretic quantities
and their estimation.

### An example: mutual information

Quantities such as "mutual information" have multiple definitions in the
literature. For example, the Shannon mutual information (MI) has both a discrete and
continuous version. For both cases, there are multiple equivalent mathematical formulas
for it. As an example, the [`MIShannon`](@ref) measure currently accepts
two definitions:

- The [`MIDefinitionGeneric`](@ref) definition computes either discrete or continuous
    MI by a sum of three entropy terms.
- The [`MIDefinitionShannonDoubleSum`](@ref) definition computes discrete Shannon MI from
    a double sum (for which joint probabilities are estimated using a contingency table).
- Even more equivalent definitions are possible, for example by directly computing some form
    of Kullback-Leibler divergence between joint distributions and products of marginal
    distributions. Each such method comes with its own advantages/disadvantages when it
    comes to practical estimation.

But this is not the only type of mutual information! Taking as a starting point some
generalized entropy definition like [`Tsallis`](@ref), several authors have proposed
variants of "Tsallis mutual information". Like Shannon MI, the Tsallis variant also has
many definitions in the scientific literature, and *not all of them are equivalent,
even though they are referred to by the same name*!
Naming ambiguities like these are likely to cause confusion.

## Designing an information measure API

To alleviate potential confusion, we have designed this package around "measures"
([`InformationMeasure`](@ref)s). Each measure is instantiated with some *definition*,
which controls which quantity is computed. Measures are combined with *estimators* that
control *how* a quantity is computed according to its definition, given som input data.

For example, the [`MIShannon`](@ref) measure can be combined with the
[`MIDefinitionGeneric`](@ref) definition, which is compatible with any
of our [`DifferentialEntropyEstimator`](@ref)s and our [`ProbabilitiesEstimator`](@ref)s.
Estimating MI from data is then simply a matter of feeding an instance of
[`MIShannon`](@ref) to [`mutualinfo`](@ref), together with an estimator and your data:

```@example
using CausalityTools
x, y = rand(10000), rand(10000)
# Use the H3-estimation method with a discrete visitation frequency based 
# probabilities estimator over a fixed grid covering the range of the data,
# which is on [0, 1].
measure = MIShannon(base = 2, definition = MIDefinitionGeneric())
est = ValueHistogram(FixedRectangularBinning(0, 1, 5))
mutualinfo(measure, est, x, y)
```

### Universal estimation

The definition-estimator based design allows us to use Julia's multiple dispatch system
to have a *single* function for computing related, but disparate concepts, like
[`mutualinfo`](@ref). Integration with probabilities estimators and entropy estimators from
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl)
is super simple, and makes it easy to develop and maintain the codebase.

With this modular API, one could estimate *any* information measure using *any* estimator.
Although the current interface doesn't allow *every* combination of measure
and estimator, you can already do a lot! To compute your favorite measure, simply find a
suitable estimator in one of the overview tables, and apply it to some input data! Follow
one of the plentiful examples for inspiration.

If you're interested in a deeper understanding, we try to give mathematical formulas and
implementation details as best we can in the docstrings of the various measures and
definitions.

### Summary

This package has been and is under heavy development. Don't hesitate to submit an
issue if you find something that doesn't work or doesn't make sense, or if there's
some functionality that you're missing!

Pull requests are also very welcome!
