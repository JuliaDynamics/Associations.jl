# CausalityTools.jl

`CausalityTools` is a Julia package that provides algorithms for *detecting 
dynamical influences* and *causal inference* based on time series data. 

For an updated overview of the field of causal inference, see for example Runge et al. (2019)[^Runge2019].

## Getting started

Examples showcasing how to use the different methods work are provided under the respective documentation packages.
Currently, this package provides two types of causal inference methods: those based on *entropy estimation*, and 
those based on geometric measures.  

### Geometric methods

Geometrically based methods rely on delay reconstructions of the time series, and numerical 
properties of these delay reconstructions, to infer causal/dynamical relationships. They take as input the
time series, and the embedding parameters are given as keyword arguments.

- Cross mapping ([`crossmap`](@ref), [`convergentcrossmap`](@ref))
- Joint distance distribution ([`jdd`](@ref))
- S-measure ([`s_measure`](@ref))

### Information theoretic methods

*Entropy* based methods for causal inference take as inputs the time series in question, and an entropy estimator of choice. Additional parameters are given as keyword arguments. 

- Generalized entropy ([`genentropy`](@ref)) and fast histograms ([`binhist`](@ref))
- Mutual information ([`mutualinfo`](@ref))
- Transfer entropy ([`transferentropy`](@ref))
- Predictive asymmetry ([`predictive_asymmetry`](@ref))

We use entropy estimators from [Entropies.jl](https://github.com/JuliaDynamics/Entropies.jl). 
Which estimator should you use? See the [list of estimators](@ref estimators). 
A good choice is to start with a [`VisitationFrequency`](@ref) estimator.

## Contributions/issues

Do you want additional methods or example systems to be implemented? Make a PR to the 
master branch in the 
[CausalityTools repo](https://github.com/JuliaDynamics/CausalityTools.jl), or open an 
issue describing the desired feature.

[^Runge2019]: Runge, J., Bathiany, S., Bollt, E., Camps-Valls, G., Coumou, D., Deyle, E., Glymour, C., Kretschmer, M., Mahecha, M. D., Muñoz-Marí, J., van Nes, E. H., Peters, J., Quax, R., Reichstein, M., Scheffer, M., Schölkopf, B., Spirtes, P., Sugihara, G., Sun, J., … Zscheischler, J. (2019). Inferring causation from time series in Earth system sciences. Nature Communications, 10(1), 1–13. https://doi.org/10.1038/s41467-019-10105-3
