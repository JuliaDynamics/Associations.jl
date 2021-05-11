# CausalityTools.jl

## Measures

`CausalityTools` is a Julia package providing algorithms for detecting dynamical relations 
between variables of complex systems based on time series data.

### Information theoretic measures

- Generalized entropy ([`genentropy`](@ref)) and fast histograms ([`binhist`](@ref))
- Mutual information ([`mutualinfo`](@ref))
- Transfer entropy ([`transferentropy`](@ref))
- Predictive asymmetry ([`predictive_asymmetry`](@ref))

See also the [list of estimators](@ref estimators) for the information theoretic measures.

### Geometric measures

- Cross mapping ([`crossmap`](@ref), [`convergentcrossmap`](@ref))
- Joint distance distribution ([`jdd`](@ref))
- S-measure ([`s_measure`](@ref))

## Contributions/issues

Do you want additional methods or example systems to be implemented? Make a PR to the 
master branch in the 
[CausalityTools repo](https://github.com/JuliaDynamics/CausalityTools.jl), or open an 
issue describing the desired feature.
