# [Discretization tutorial](@id discretization_tutorial)

There are two main ways of discretizing data in CausalityTools. They are implemented as 
the [`CodifyPoints`](@ref) and [`CodifyVariables`](@ref) types, which are used as 
input to the [`codify`](@ref) function (extended from ComplexityMeasures.jl to multiple 
variables).

## [Encoding *rows* (one *point* at a time)](@id tutorial_codify_points)

In some cases, it may be desireable to encode data on a row-wise basis. This 
typically happens when working with pre-embedded time series. If we want to 
apply something like [`OrdinalPatternEncoding`](@ref) to a pre-embedded 
[`StateSpaceSet`](@ref), then we must encode each *point* individually,
respecting the fact that time ordering is already taken care of by the 
embedding procedure.

```@example
using CausalityTools
using Random; rng = Xoshiro(1234)

# The first variable is 2-dimensional and has 50 points
x = StateSpaceSet(rand(rng, 50, 2))
# The second variable is 3-dimensional and has 50 points
y = StateSpaceSet(rand(rng, 50, 3))
# The third variable is 4-dimensional and has 50 points
z = StateSpaceSet(rand(rng, 50, 4))

ox = OrdinalPatternEncoding(2)
oy = OrdinalPatternEncoding(3)
oz = OrdinalPatternEncoding(4)

# This given three column vectors of integers.
# `x` has been encoded usign `ox` on a point-by-point basis
# `y` has been encoded using `oy` on a point-by-point basis
# `z` has been encoded using `oz` on a point-by-point basis
cx, cy, cz = codify(CodifyPoints(ox, oy, oz), x, y, z)

[cx cy cz]
```

Notice that the 2-dimensional `x` has been encoded into integer values `1` or `2`, because
there are `2!` possible ordinal patterns for dimension `m = 2`. The 3-dimensional `y` has 
been encoded into integers in the range `1` to `3! = 6`, while the 4-dimensional `z` is 
encoded into an even larger range of integers, because the number of possible ordinal patterns
is `4! = 24` for 4-dimensional embedding vectors.

## Encoding *columns* (one variable at a time)

Sometimes, it may be desireable to encode input data one variable/column at a time.
This typically happens when the input are either a single or multiple timeseries.

To encode columns, we apply an [`Encoding`](@ref) using a sliding window across each input variable. The width of the window is determined by the chosen encoding.
For example, using [`ValueBinning`](@ref) will encode `N` value into `N` discretized
values.

```@example
using CausalityTools
using Random; rng = Xoshiro(1234)

x = rand(rng, 100)
o = ValueBinning(3)
cx = codify(CodifyVariables(o), x)
```

We can verify that [`ValueBinning`](@ref) preserves the cardinality of the input dataset.

```@example
length(x) == length(cx)
```

Other outcome spaces such as [`Dispersion`](@ref) or [`OrdinalPatterns`](@ref) do not 
preserve the cardinality of the input dataset, because when applied in a sliding window,
they compress embedding vectors into single integers. This means that some points at the 
end of each input variable are lost.

```@example
using CausalityTools
using Random; rng = Xoshiro(1234)

x = rand(rng, 100)
o = OrdinalPatterns(m = 3)
cx = codify(CodifyVariables(o), x)
```

We can also simultaneously encode each variable/column of a [`StateSpaceSet`](@ref), as long 
as we apply an encoding that results in the *same* number of encoded data points.

```@example
using CausalityTools
using Random; rng = Xoshiro(1234)

x = rand(rng, 100)
y = rand(rng, 100)
o = OrdinalPatterns(m = 3)
# Alternatively provide a tuple of input time series: codify(CodifyVariables(o), (x, y))
cx, cy = codify(CodifyVariables(o), StateSpaceSet(x, y)) 

[cx cy]
```