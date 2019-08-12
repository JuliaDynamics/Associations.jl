# Transfer entropy estimation

For high-level methods for TE estimation, see 
[# Convenience functions for TE estimation](@ref). Some methods for estimating TE from 
two time series, potentially conditioned on a third time series, are provided. The methods 
below offers completely control over the estimation procedure.

## Rectangular partitions

The general workflow for estimating transfer entropy over rectangular partitions is as follows.

```@docs
transferentropy(pts, vars::TEVars, ϵ::RectangularBinning, estimator::TransferEntropyEstimator; b = 2)
```

Valid estimator types are `VisitationFrequency` and `TransferOperatorGrid`, so to estimate 
transfer entropy, you would use the following methods.

```@docs
transferentropy(pts, vars::TEVars, ϵ::RectangularBinning, estimator::VisitationFrequency; b = 2)
```

```@docs
transferentropy(pts, vars::TEVars, ϵ::RectangularBinning, estimator::TransferOperatorGrid; b = 2)
```

## Triangulated partitions

```@docs
transferentropy(μ::AbstractTriangulationInvariantMeasure, vars::TEVars,
        binning_scheme::RectangularBinning;
        estimator = VisitationFrequency(), n::Int = 10000, b = 2)
```