# Transfer entropy (TE) estimators

## Common interface

We provide a common `transferentropy(driver, response; kwargs...)` function,
which may be nice for initial exploration of your data.

*Note: the common interface uses default values for estimator parameters that
may not be suited for your data. For real applications, it is highly
recommended to use the underlying estimators directly.*

## Documentation

```@docs
te
```
