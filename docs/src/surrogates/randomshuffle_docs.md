# Surrogate methods

## Random shuffle surrogates 


### Valid inputs
Random shuffle surrogates may be generated from the following inputs:

- `AbstractArray{T, 1}` instances (scalar-valued data series)
- `AbstractArray{Number, 2}` instances (multivarate scalar-valued data series), for which surrogates are generated column-wise. 
- `Dataset` instances from [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl), for which surrogates are generated column-wise. 
- `Embedding` instances, for which surrogates are generated variable-wise (row-wise on the points).


## Documentation

```@docs
randomshuffle
```