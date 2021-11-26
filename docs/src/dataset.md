
# Multivariate [`Dataset`](@ref)s

Univariate timeseries are given as `AbstractVector{<:Real}`. Multivariate time series can
be represented by the `Dataset`-type from [`DelayEmbeddings.jl`](https://github.com/JuliaDynamics/DelayEmbeddings.jl), where each observation is a D-dimensional data point represented by a static vector. See the [`DynamicalSystems.jl` documentation](https://juliadynamics.github.io/DynamicalSystems.jl/dev/) for more info.

```@docs
Dataset
```
