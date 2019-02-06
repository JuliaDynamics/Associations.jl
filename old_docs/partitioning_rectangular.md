# Rectangular partitions
Currently, there are four different ways of constructing rectangular partitions.
Functions using partitions take the `ϵ` argument, and the following `ϵ` will
work:

* `ϵ::Int` divide each axis into `ϵ` intervals of the same size.
* `ϵ::Float` divide each axis into intervals of size `ϵ`.
* `ϵ::Vector{Int}` divide the i-th axis into `ϵᵢ` intervals of the same size.
* `ϵ::Vector{Float64}` divide the i-th axis into intervals of size `ϵᵢ`.

Below, we plot the box coverings resulting from each of the methods on a
set of random points (here, we represent the points as an embedding, but
creating box coverings of a set of points works too).

```@setup rectangular_partitions
using CausalityTools
using Plots
```

```@example rectangular_partitions
E = customembed(rand(3, 100))
```

## Hyper-rectangles by subdivision of axes (`ϵ::Int`)

Subdivide each axis into ten intervals of equal length:

```@example rectangular_partitions
```
