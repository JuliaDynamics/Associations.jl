# Bivariate and multivariate information measures

Information measures are functionals of probability mass functions (PMFs),
or probability densities (following the convention in
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl)).

While ComplexityMeasures.jl implements univariate measures, we here focus on bivariate
and multivariate measures. Analogously to the univariate case, these are simply
functionals of multidimensional  probability distributions/densities.

## Encoding API

```@docs
encoding(::PointEncoding{N}, x::Vararg{<:AbstractStateSpaceSet, N})
```
