# Quickstart TE estimators

```@setup s
using CausalityTools
# Generate two time series and embed using E = {(x2(t+1), x2(t), x1(t))}
ts_length = 10
ts1, ts2 = ([diff(rand(ts_length)) for i = 1:2]...)
E = embed([ts1, ts2], [2, 2, 1], [1, 0, 0])
ϵ = 2

# Which variables of the embedding goes into which marginals?
vars  = TEVars([1], [2], [3], Int[])

# Compute the transfer entropy
transferentropy_visitfreq(E, ϵ, vars)
```

## Visitation frequency based estimator (traditional TE)

Embed two time series, measure the transfer entropy from ``x1`` to ``x2`` using a partition of scale `ϵ = 3` in dimension ``3`` (giving ``3^3`` equidistantly spaced rectangular bins).

```@example s
# Generate two uncoupled red noise processes time series
ts_length = 100
ts1, ts2 = ([cumsum(rand(ts_length)) for i = 1:2]...)

# Embed using E = {(ts2(t+1), ts2(t), ts1(t))}
E = embed([ts1, ts2], [2, 2, 1], [1, 0, 0])

# The partition scheme
ϵ = 3

# Which variables of the embedding goes into which marginals?
vars  = TEVars([1], [2], [3], Int[])

# Compute the transfer entropy
tefreq(E, ϵ, vars)
```

## Transfer operator visitation frequency based estimator

Embed two time series, measure the transfer entropy from ``x1`` to ``x2`` using a partition of scale `ϵ = 3` in dimension ``3`` (giving ``3^3`` equidistantly spaced rectangular bins).

```@example s
# Generate two uncoupled red noise processes time series
ts_length = 100
ts1, ts2 = ([cumsum(rand(ts_length)) for i = 1:2]...)

# Embed using E = {(ts2(t+1), ts2(t), ts1(t))}
E = embed([ts1, ts2], [2, 2, 1], [1, 0, 0])

# The partition scheme
ϵ = 3

# Which variables of the embedding goes into which marginals?
vars  = TEVars([1], [2], [3], Int[])

# Compute the transfer entropy
tetofreq(E, ϵ, vars)
```
