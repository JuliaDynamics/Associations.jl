# Generalized embedding and data structures

Most of the causal inference methods in CausalityTools rely internally on the concept of 
generalized delay reconstructions.  

To construct these, we use `GeneralizedEmbedding` and `genembed` from [`DelayEmbeddings.jl`](https://github.com/JuliaDynamics/DelayEmbeddings.jl). Delay reconstructions are stored internally as `Dataset`s.

```@docs 
GeneralizedEmbedding
genembed
Dataset
```