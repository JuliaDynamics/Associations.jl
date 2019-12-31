# [Causal analysis](@id causal_analysis)

A [`CausalAnalysis`](@ref) is a wrapper type that stores the result of a call to the 
`causality` method. There are two ways of calling `causality`: 

- Input data are provided as the first arguments, and the causality test as the last 
    argument. This will immediately return the raw values for the causality statistic.
- The causality test is provided as the first argument, and the input data after. In
    this case, the test parameters, the input data and the values for the causality 
    statistic computed on the data will be wrapper in a [`CausalAnalysis`](@ref) 
    instance.

```@docs
CausalAnalysis
```

```@docs
summarise
```
