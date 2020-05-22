# Estimating transfer entropy directly from time series

The following methods allows you to skip the technicalities of constructing
generalised delay reconstructions and taking care of the mapping of marginals.

## [Transfer entropy between two time series](@id te_timeseries)

```@docs
transferentropy(source::AbstractArray{<:Real, 1},
        response::AbstractArray{<:Real, 1},
        k::Int, l::Int, m::Int,
        binning_scheme::Union{RectangularBinning, Vector{RectangularBinning}})
```

## [Transfer entropy between two time series conditioned on a third time series](@id te_timeseries_cond)

```@docs
transferentropy(source::AbstractArray{<:Real, 1},
        response::AbstractArray{<:Real, 1},
        cond::AbstractArray{<:Real, 1},
        k::Int, l::Int, m::Int, n::Int,
        binning_scheme::Union{RectangularBinning, Vector{RectangularBinning}})
```
