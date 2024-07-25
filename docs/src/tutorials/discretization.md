# [Discretization tutorial](@id extended_example_discretization)

When discretizing, what happens is that we "encode" input data into an intermediate representation indexed by the positive integers. This intermediate representation is called an "encoding". This is useful in several ways:

- Once a dataset has been encoded into integers, we can estimate [`Counts`](@ref) or [`Probabilities`](@ref) ([tutorial](@ref tutorial_probabilities)).
- Once probabilities have been estimated, one can use these to estimate [`MultivariateInformationMeasure`](@ref) ([tutorial](@ref tutorial_infomeasures)).

## Two ways of encoding input data

There are two main ways of discretizing data in CausalityTools. 
- The [`CodifyPoints`](@ref) discretization scheme encodes input data on a point-by-point 
    basis by applying some [`Encoding`](@ref) to each point.
- The [`CodifyVariables`](@ref) discretization scheme encodes input data on a column-by-column
    basis by applying a sliding window to each column, and encoding the data within the sliding window according to some [`OutcomeSpace`](@ref) (*Internally, this uses [`codify`](@ref)*).

!!! note 
    [`Encoding`](@ref), [`OutcomeSpace`](@ref) and [`codify`](@ref) are all from
    [ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl),
    but is here used to discretize multiple input variables instead of just one input
    variable.

### [Encoding *rows* (one *point* at a time)](@id tutorial_codify_points)

In some cases, it may be desireable to encode data on a row-wise basis. This 
typically happens when working with pre-embedded time series. If we want to 
apply something like [`OrdinalPatternEncoding`](@ref) to a pre-embedded 
[`StateSpaceSets.StateStateSpaceSet`](@extref), then we must encode each *point* individually,
respecting the fact that time ordering is already taken care of by the 
embedding procedure. [`CodifyPoints`](@ref) ensures input data are encoded 
on a point-by-point basis.

```@example example_encode_points
using CausalityTools
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

### Encoding *columns* (one variable at a time)

Sometimes, it may be desireable to encode input data one variable/column at a time.
This typically happens when the input are either a single or multiple timeseries.

To encode columns, we apply an [`Encoding`](@ref) using a sliding window across each input variable. 
The width of the window is determined by the chosen encoding.
For example, using [`ValueBinning`](@ref) will encode `N` value into `N` discretized
values. [`CodifyVariables`](@ref) is used to enforce a sliding window encoding on a 
per-variable basis.

```@example example_encode_vars
using CausalityTools
using Random; rng = Xoshiro(1234)

x = rand(rng, 100)
o = ValueBinning(3)
cx = codify(CodifyVariables(o), x)
```

We can verify that [`ValueBinning`](@ref) preserves the cardinality of the input dataset.

```@example example_encode_vars
length(x) == length(cx)
```

Other outcome spaces such as [`Dispersion`](@ref) or [`OrdinalPatterns`](@ref) do not 
preserve the cardinality of the input dataset, because when applied in a sliding window,
they compress sliding windows consisting of potentially multiple points into single integers. This means that some points at the 
end of each input variable are lost.

```@example example_encode_vars
using CausalityTools
using Random; rng = Xoshiro(1234)

x = rand(rng, 100)
o = OrdinalPatterns(m = 3)
cx = codify(CodifyVariables(o), x)
```

We can also simultaneously encode each variable/column of a [`StateSpaceSet`](@ref), as long 
as we apply an encoding that results in the *same* number of encoded data points.

```@example example_encode_vars
using CausalityTools
using Random; rng = Xoshiro(1234)

x = rand(rng, 100)
y = rand(rng, 100)
o = OrdinalPatterns(m = 3)
# Alternatively provide a tuple of input time series: codify(CodifyVariables(o), (x, y))
cx, cy = codify(CodifyVariables(o), StateSpaceSet(x, y)) 

[cx cy]
```

