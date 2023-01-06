# Design philosophy

## Fundamentals

The building blocks of most information theoretic measures are [probabilities](@ref probabilities_header) and [entropies](@ref entropy_header). However, these quantities
also appear in many other applications, for example in complexity quantification.

Realizing this, we decided the gather the fundamentals required for estimation of
probabilities and entropies in a single package:
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl).
ComplexityMeasures.jl was built with modularity in mind, and provides a plethora of
estimators of probabilities and generalized entropies, both discrete and continuous.
These estimators are used frequently throughout CausalityTools.jl.

## Definitions: the same name for different things

In contrast to generalized entropies, which each have *one* definition, it gets a bit more
complicated when it comes to the "higher-level" measures we provide here.

Upon doing a literature review on the possible variants of information theoretic measures,
it become painstakingly obvious that authors use *the same name for different concepts*.
For novices in the field of information theory and its practical applications, this
can be very confusing. This package is designed to alleviate any confusion regarding
the names of information theoretic quantities and their estimation.

### Example: mutual information definitions

Quantities such as "mutual information" have multiple definitions in the
literature. For example, the Shannon mutual information (MI; [`MIShannon`](@ref)) has
both a discrete and continuous version, and there there are multiple equivalent mathematical formulas for them. Thus, there are also many ways to estimate these quantities.

But this is not the only type of mutual information! Taking as a starting point some
generalized entropy definition like [`Tsallis`](@ref), several authors have proposed
variants of "Tsallis mutual information". Like Shannon MI, the Tsallis variant also has
many definitions in the scientific literature, and *not all of them are equivalent,
even though they are referred to by the same name*! Naming ambiguities like these are
likely to cause confusion.

To alleviate potential confusion, we name information theoretic measures in this
package starting with a few letters indicating which *type* of measure it is, for example
`MI` for mutual information and `CE` for conditional entropy. Onto that we append
a definition name, like `Shannon`, `Tsallis` or `Renyi`. Finally, if there are multiple
non-equivalent definitions, we append an author name to the measure to specify 
which definition is to be computed. For example, [`MIShannon`](@ref) only exists
as one measure, because all of its definitions are equivalent. But
[`MITsallisFuruichi`](@ref) and [`MITsallisMartin`](@ref) are separate measures, 
because they represent different mathematical formulas.

Measures are combined with *estimators* that control *how* a quantity is computed
according to its definition, given som input data.

For example, to [`MIShannon`](@ref) can be combined with
any [`DifferentialEntropyEstimator`](@ref) to compute differential Shannon mutual
information. For the discrete quantity, you can combine [`MIShannon`](@ref) with
any compatible [`ProbabilitiesEstimator`](@ref)s.

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
# Use the H3-estimation method with a discrete visitation frequency based 
# probabilities estimator over a fixed grid covering the range of the data,
# which is on [0, 1].
measure = MIShannon(base = 2)
est = ValueHistogram(FixedRectangularBinning(0, 1, 5))
mutualinfo(measure, est, x, y)
```

This is in fact is equivalent to using a [`ContingencyMatrix`](@ref):

```@example mi_demonstration
c = contingency_matrix(est, x, y)
mutualinfo(measure, c)
```

However, the [`ContingencyMatrix`](@ref) can also be used with categorical data.

### Universal estimation

The measure-estimator based design allows us to use Julia's multiple dispatch system
to have a *single* function for computing a single concept in multiple ways.
Integration with probabilities estimators and entropy estimators from
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl)
is super simple, and makes it easy to develop and maintain the codebase.

With this modular API, one could in principle estimate *any* information measure using *any* estimator. Although the current interface doesn't allow *every* combination of measure
and estimator, and it's probably not theoretically meaningful,
you can already do a lot! To compute your favorite measure, simply find a
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
