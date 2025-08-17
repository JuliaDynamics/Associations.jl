
```@meta
CollapsedDocStrings = true
```

# [Multivariate counts and probabilities API](@id counts_and_probabilities_api)

For counting and probabilities, Associations.jl extends the single-variable machinery
in ComplexityMeasures.jl to multiple variables.

```@docs
Associations.Counts
Associations.counts(::OutcomeSpace)
```

```@docs
Associations.Probabilities
Associations.probabilities(::OutcomeSpace)
```

The utility function [`marginal`](@ref) is also useful.

```@docs
marginal
```

## [Example: estimating [`Counts`](@extref ComplexityMeasures.Counts) and [`Probabilities`](@extref ComplexityMeasures.Probabilities)](@id tutorial_probabilities)

Estimating multivariate counts (contingency matrices) and PMFs is simple. If the data are pre-discretized, then
we can use [`UniqueElements`](@extref ComplexityMeasures.UniqueElements) to simply count the number of occurrences.

```@example counts_probs_tutorial
using Associations
n = 50 # the number of samples must be the same for each input variable
x = rand(["dog", "cat", "snake"], n)
y = rand(1:4, n)
z = rand([(2, 1), (0, 0), (1, 1)], n)
discretization = CodifyVariables(UniqueElements())
counts(discretization, x, y, z)
```

Probabilities are computed analogously, except counts are normalized to sum to `1`.

```@example counts_probs_tutorial
discretization = CodifyVariables(UniqueElements())
probabilities(discretization, x, y, z)
```

For numerical data, we can estimate both counts and probabilities using [`CodifyVariables`](@ref)
with any count-based [`OutcomeSpace`](@extref ComplexityMeasures.OutcomeSpace).

```@example counts_probs_tutorial
using Associations
x, y = rand(100), rand(100)
discretization = CodifyVariables(BubbleSortSwaps(m = 4))
probabilities(discretization, x, y)
```

For more fine-grained control, we can use [`CodifyPoints`](@ref) with one or several [`Encoding`](@extref ComplexityMeasures.Encoding)s.

```@example counts_probs_tutorial
using Associations
x, y = StateSpaceSet(rand(1000, 2)), StateSpaceSet(rand(1000, 3))

 # min/max of the `rand` call is 0 and 1
precise = true # precise bin edges
r = range(0, 1; length = 3)
binning = FixedRectangularBinning(r, dimension(x), precise)
encoding_x = RectangularBinEncoding(binning, x)
encoding_y = CombinationEncoding(RelativeMeanEncoding(0.0, 1, n = 2), OrdinalPatternEncoding(3))
discretization = CodifyPoints(encoding_x, encoding_y)

# now estimate probabilities
probabilities(discretization, x, y)
```
