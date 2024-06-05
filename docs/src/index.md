![CausalityTools.jl static logo](assets/logo-large.png)

```@docs
CausalityTools
```

## Goals

Causal inference, and quantification of association in general, is fundamental to
most scientific disciplines. There exists a multitude of bivariate and multivariate
association measures in the scientific literature. However, beyond the most basic measures,
most methods aren't readily available for practical use. Most scientific papers don't
provide code, which makes reproducing them difficult or impossible, without
investing significant time and resources into deciphering and understanding the original
papers to the point where an implementation is possible. To make reliable inferences,
proper independence tests are also crucial.

Our main goal with this package is to provide an easily extendible library of
association measures, a as-complete-as-possible set of their estimators.

On top of the association measures, we've built a framework for independence testing. 
Finally, at the highest level, we provide a framework for independence testing, 
where you can freely combine any compatible association measure and independence test 
to deduce causal graphs.

This package also extends the univariate information measure API in
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl).
