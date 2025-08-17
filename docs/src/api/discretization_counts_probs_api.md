```@meta
CollapsedDocStrings = false
```

# Discretization API

## Encoding multiple input datasets

A fundamental operation when computing multivariate information measures from data is *discretization*. 
When discretizing, what happens is that we "encode" input data into an intermediate representation indexed by the positive integers. This intermediate representation is called an "encoding". This is useful in several ways:

- Once a dataset has been encoded into integers, we can estimate [`Counts`](@extref ComplexityMeasures.Counts) or [`Probabilities`](@extref ComplexityMeasures.Probabilities) ([tutorial](@ref tutorial_probabilities)).
- Once probabilities have been estimated, one can use these to estimate [`MultivariateInformationMeasure`](@ref) ([tutorial](@ref tutorial_infomeasures)).

 The following functions and types are used by Associations.jl to perform discretization of input data.

```@docs
Discretization
CodifyVariables
CodifyPoints
codify
```

In summary, the two main ways of discretizing data in Associations are as follows.

- The [`CodifyPoints`](@ref) discretization scheme encodes input data on a point-by-point 
    basis by applying some [`Encoding`](@extref ComplexityMeasures.Encoding) to each point.
- The [`CodifyVariables`](@ref) discretization scheme encodes input data on a column-by-column
    basis by applying a sliding window to each column, and encoding the data within the sliding window according to some [`OutcomeSpace`](@extref ComplexityMeasures.OutcomeSpace) (*Internally, this uses [`codify`](@ref)*).

!!! note 
    [`Encoding`](@extref ComplexityMeasures.Encoding), [`OutcomeSpace`](@extref ComplexityMeasures.OutcomeSpace) and [`codify`](@ref) are all from
    [ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl).
    In this package, they are used to discretize multiple input variables instead of just one input
    variable.


### Encoding per point/row

In some cases, it may be desireable to encode data on a row-wise basis. This 
typically happens when working with pre-embedded time series or [`StateSpaceSet`](@extref StateSpaceSets.StateSpaceSet)s 
(respecting the fact that time ordering is already taken care of by the 
embedding procedure). 
If we want to apply something like [`OrdinalPatternEncoding`](@extref ComplexityMeasures.OrdinalPatternEncoding) to such data, then 
we must encode each *point* individually, such that vectors like `[1.2, 2.4, 4.5]` or 
`["howdy", "partner"]` gets mapped to an integer. The [`CodifyPoints`](@ref) discretization 
intstruction ensures input data are encoded on a point-by-point basis.

A point-by-point discretization using [`CodifyPoints`](@ref) is formally done by applying some [`Encoding`](@extref ComplexityMeasures.Encoding) to each input data point. You can pick between the following encodings, or combine 
them in arbitrary ways using [`CombinationEncoding`](@extref ComplexityMeasures.CombinationEncoding).

- [`Encoding`](@extref ComplexityMeasures.Encoding)
- [`GaussianCDFEncoding`](@extref ComplexityMeasures.GaussianCDFEncoding)
- [`OrdinalPatternEncoding`](@extref ComplexityMeasures.OrdinalPatternEncoding)
- [`RelativeMeanEncoding`](@extref ComplexityMeasures.RelativeMeanEncoding)
- [`RelativeFirstDifferenceEncoding`](@extref ComplexityMeasures.RelativeFirstDifferenceEncoding)
- [`UniqueElementsEncoding`](@extref ComplexityMeasures.UniqueElementsEncoding)
- [`RectangularBinEncoding`](@extref ComplexityMeasures.RectangularBinEncoding)
- [`CombinationEncoding`](@extref ComplexityMeasures.CombinationEncoding)

#### [Examples: encoding *rows* (one *point* at a time)](@id tutorial_codify_points)

We'll here use the [`OrdinalPatternEncoding`](@extref ComplexityMeasures.OrdinalPatternEncoding) with differing parameter `m` to encode 
multiple [`StateSpaceSet`](@extref StateSpaceSets.StateSpaceSet) of differing dimensions.

```@example example_encode_points
using Associations
using StateSpaceSets
using Random; rng = Xoshiro(1234)

# The first variable is 2-dimensional and has 50 points
x = StateSpaceSet(rand(rng, 50, 2))
# The second variable is 3-dimensional and has 50 points
y = StateSpaceSet(rand(rng, 50, 3))
# The third variable is 4-dimensional and has 50 points
z = StateSpaceSet(rand(rng, 50, 4))

# One encoding scheme per input variable
# encode `x` using `ox` on a point-by-point basis (Vector{SVector{4}} → Vector{Int})
# encode `y` using `oy` on a point-by-point basis (Vector{SVector{3}} → Vector{Int})
# encode `z` using `oz` on a point-by-point basis (Vector{SVector{2}} → Vector{Int})
ox = OrdinalPatternEncoding(2)
oy = OrdinalPatternEncoding(3)
oz = OrdinalPatternEncoding(4)

# This given three column vectors of integers.
cx, cy, cz = codify(CodifyPoints(ox, oy, oz), x, y, z)

[cx cy cz]
```

Notice that the 2-dimensional `x` has been encoded into integer values `1` or `2`, because
there are `2!` possible ordinal patterns for dimension `m = 2`. The 3-dimensional `y` has 
been encoded into integers in the range `1` to `3! = 6`, while the 4-dimensional `z` is 
encoded into an even larger range of integers, because the number of possible ordinal patterns
is `4! = 24` for 4-dimensional embedding vectors.

### Encoding per variable/column

Sometimes, it may be desireable to encode input data one variable/column at a time.
This typically happens when the input are either a single or multiple timeseries.

To encode columns, we move a sliding window across each input variable/column and 
encode points within that window. Formally, such a sliding-window discretization 
is done by using the [`CodifyVariables`](@ref) discretization scheme, which takes
as input some [`OutcomeSpace`](@extref ComplexityMeasures.OutcomeSpace) that dictates how each window is encoded, and 
also dictates the width of the encoding windows. 

For column/variable-wise encoding, you can pick between the following outcome spaces.
- [`OutcomeSpace`](@extref ComplexityMeasures.OutcomeSpace)
- [`UniqueElements`](@extref ComplexityMeasures.UniqueElements)
- [`CosineSimilarityBinning`](@extref ComplexityMeasures.CosineSimilarityBinning)
- [`Dispersion`](@extref ComplexityMeasures.Dispersion)
- [`OrdinalPatterns`](@extref ComplexityMeasures.OrdinalPatterns)
- [`OrdinalPatterns`](@extref ComplexityMeasures.OrdinalPatterns)
- [`BubbleSortSwaps`](@extref ComplexityMeasures.BubbleSortSwaps)
- [`ValueBinning`](@extref ComplexityMeasures.ValueBinning)
- [`RectangularBinning`](@extref ComplexityMeasures.RectangularBinning)
- [`FixedRectangularBinning`](@extref ComplexityMeasures.FixedRectangularBinning)


#### Example: encoding *columns* (one variable at a time)

Some [`OutcomeSpace`](@extref ComplexityMeasures.OutcomeSpace)s dictate a sliding window which has the width of one element
when used with [`CodifyVariables`](@ref). [`ValueBinning`](@extref ComplexityMeasures.ValueBinning) is such an outcome space.

```@example example_encode_vars
using Associations
using Random; rng = Xoshiro(1234)

x = rand(rng, 100)
o = ValueBinning(3)
cx = codify(CodifyVariables(o), x)
```

We can verify that [`ValueBinning`](@extref ComplexityMeasures.ValueBinning) preserves the cardinality of the input dataset.

```@example example_encode_vars
length(x) == length(cx)
```

Other outcome spaces such as [`Dispersion`](@extref ComplexityMeasures.Dispersion) or [`OrdinalPatterns`](@extref ComplexityMeasures.OrdinalPatterns) do not 
preserve the cardinality of the input dataset when used with [`CodifyVariables`](@ref). This is 
because when they are applied in a sliding window, they compress sliding windows consisting of 
potentially multiple points into single integers. This means that some points at the 
end of each input variable are lost. For example, with [`OrdinalPatterns`](@extref ComplexityMeasures.OrdinalPatterns), the number 
of encoded points decrease with the embedding parameter `m`.

```@example example_encode_vars
using Associations
using Random; rng = Xoshiro(1234)

x = rand(rng, 100)
o = OrdinalPatterns(m = 3)
cx = codify(CodifyVariables(o), x)
```

We can simultaneously encode multiple variable/columns of a [`StateSpaceSet`](@extref StateSpaceSets.StateSpaceSet) using 
the same outcome space, as long as the operation will result in the *same* number of encoded 
data points for each column.

```@example example_encode_vars
using Associations
using Random; rng = Xoshiro(1234)

x = rand(rng, 100)
y = rand(rng, 100)
o = OrdinalPatterns(m = 3)
# Alternatively provide a tuple of input time series: codify(CodifyVariables(o), (x, y))
cx, cy = codify(CodifyVariables(o), StateSpaceSet(x, y)) 

[cx cy]
```

