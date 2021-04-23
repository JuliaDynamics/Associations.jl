# CausalityTools.jl

`CausalityTools` is a Julia package providing algorithms for detecting dynamical relations 
between variables of complex systems based on time series data.

## Causal inference methods for time series

### Information theoretic measures

- Generalized entropy ([`genentropy`](@ref)) and fast histograms ([`binhist`](@ref))
- Mutual information ([`mutual_info`](@ref))
- Transfer entropy ([`transferentropy`](@ref))
- Predictive asymmetry ([`predictive_asymmetry`](@ref))

See also the [list of estimators](@ref estimators) for the information theoretic measures.

### Geometric measures

- Cross mapping ([`crossmap`](@ref), [`convergentcrossmap`](@ref))
- Joint distance distribution ([`jdd`](@ref))
- S-measure ([`s_measure`](@ref))

## Examples of coupled dynamical systems

- [Discrete systems](@ref discrete_systems)
- [Continuous systems](@ref continuous_systems)

## Contributions

Do you want additional methods or example systems to be implemented? Make a PR to the 
master branch in the 
[CausalityTools repo](https://github.com/JuliaDynamics/CausalityTools.jl), or open an 
issue describing the desired feature.
