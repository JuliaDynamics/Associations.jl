# PerronFrobenius.jl

## Transfer operator approximation 

The `transferoperator` function is used for approximating the transfer operator 
over a partition of the input data. `transferoperator` dispatches on subtypes of 
`TransferOperator`, which controls the partition and the method used for estimating 
transition probabilities.

```@docs
transferoperator
```

```@docs 
SingleGrid
SimplexPoint
SimplexExact
```

## Invariant measure approximation

From the transfer operators, one may compute invariant probability distributions over the 
partitions too. `invariantmeasure` also dispatches on subtypes of 
`TransferOperator`, but can also be used on pre-computed from pre-comptued transfer operator 
approximations.

```@docs
invariantmeasure
```
