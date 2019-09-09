
# [Custom delay reconstruction](@id custom_delay_reconstruction)

To create delay embeddings from sequential data, use `customembed`. This function returns allows for more flexibility than `embed` in DynamicalSystems.jl[^1], but also returns a `Dataset`, so these embedding methods may be used interchangeably.

## Documentation
```@docs 
customembed
```

```@docs 
Positions
```

```@docs 
Lags
```

[^1]:
    Datseris, George. "DynamicalSystems. jl: A Julia software library for chaos and nonlinear dynamics." J. Open Source Software 3.23 (2018): 598.