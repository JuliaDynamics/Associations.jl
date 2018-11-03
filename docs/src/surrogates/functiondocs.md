# Function overview

```@setup s
using CausalityTools
using Plots
```

## Constrained realizations

The following surrogate functions generates constrained realizations of their
input data series `ts`, i.e. values of `ts` are shuffled in a way that preserves
or destroys certain statistical properties of `ts`.

```@docs
randomshuffle(ts)
```

```@docs
aaft(ts)
```

```@docs
iaaft(ts)
```

## Unconstrained realizations

```@docs
randomphases(ts)
```

```@docs
randomamplitudes(ts)
```
