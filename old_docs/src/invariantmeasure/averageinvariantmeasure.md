# Invariant measure estimation

## Rectangular grid, averaging over multiple parition resolutions

`AverageRectangularInvariantMeasure` acts both as a type holding all information used to 
estimate an invariant measure for a rectangular partition induced by multiple other
rectangular partitions, and as a constructor. 

The constructor takes three arguments: the set of points from which to estimate the measure, 
and a target [binning scheme](../../glossary/dicretization.md) giving the final resolution `ϵF`, 
and a vector of source [binning schemes](../../glossary/dicretization.md) `ϵjs` from which we 
induce and average the measure.

### Example 

```julia 
# Create some random points
pts = [rand(3) for i = 1:20]

# Find the invariant measure using a rectangular partition with edge lengths derived from
# deviding all coordinate axes into 4 equally spaced intervals.
# However, don't consider the measure induced directly at this resolution. Rather, consider 
# the partitions formed by first dividing all axes into 5 intervals, then dividing each axis 
# into 3, 4 and 6 intervals, and finally dividing all axes into 9 intervals,
# then consider the measure induced by these finer partitions onto the coarser partition.
avg_measure = averagerectangularinvariantmeasure(pts, 6, [5, [7, 4, 6], 9])
```

There's a simple plot recipe for visualizing the visited bins at different resolutions:

```julia
# Plot the boxes with nonzero measure, along with the points the measure was estimated from
plot(avg_measure)
```

![](average_invariant_measure_4_5-346-6_20points.svg)

### Documentation

```@docs
averagerectangularinvariantmeasure
```
