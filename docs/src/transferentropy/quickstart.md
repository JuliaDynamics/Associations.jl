# Quickstart TE estimators

```@setup s
using CausalityTools
```

## Visitation frequency based estimator (traditional TE)

Embed two time series, measure the transfer entropy from ``x1`` to ``x2`` using a partition of scale `ϵ = 3` in dimension ``3`` (giving ``3^3`` equidistantly spaced rectangular bins).

```@example s
# Generate time series and embed using E = {(x2(t+1), x2(t), x1(t))}
ts_length = 30
ts1, ts2 = ([diff(rand(ts_length)) for i = 1:2]...)
E = embed([ts1, ts2], [2, 2, 1], [1, 0, 0])
ϵ = 3

# Which variables of the embedding goes into which marginals?
vars  = TEVars([1], [2], [3], Int[])

# Compute the transfer entropy
transferentropy_visitfreq(E, ϵ, vars)
```
