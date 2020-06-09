# Predictive asymmetry

CausalityTools provide the following interface for computing the predictive asymmetry (Haaga et al., 2020).
The `predictive_asymmetry` method dispatches on different [transfer entropy estimators](@ref te_estimators),
(the algorithm uses a form of transfer entropy under the hood).

```@docs
predictive_asymmetry
```