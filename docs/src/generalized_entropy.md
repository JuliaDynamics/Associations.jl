# Generalized entropy

The following two functions are used for probability and entropy estimation:

* [`probabilities`](@ref) which computes probability distributions of given datasets
* [`genentropy`](@ref) which uses the output of [`probabilities`](@ref), or a set of
    pre-computed [`Probabilities`](@ref), to calculate entropies.

See the [Entropies.jl](https://github.com/JuliaDynamics/Entropies.jl) documentation for 
details.

```@docs
Entropies.genentropy
```

```@docs
Probabilities
probabilities
probabilities!
ProbabilitiesEstimator
```

```@docs
binhist
```
