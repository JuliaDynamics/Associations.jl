# Transfer entropy (TE) and information flows

```@setup s
using CausalityTools
using DynamicalSystems

# Generate data
s = CausalityTools.Systems.logistic3()
orbit = trajectory(s, 200)
x1, x2 = orbit[:, 1], orbit[:, 2]

# Which of the time series do we embed?
which_ts = [x1, x2] # these are now indexed 1, 2

# Which variables of the embedding will each of those
# time series correspond to? And what are the lags?
which_pos = [2, 2, 1, 1]
which_lags = [1, 0, 0, -1]
E = embed(which_ts, which_pos, which_lags)

# Construct the TEVars instance.
v = TEVars([1], [2], [3, 4], Int[])

# Subdivide each axis into 15 intervals of equal size, resulting in
# rectangular boxes.
ϵ = 15

# Estimate transfer entropy at scale `ϵ`.
transferentropy_transferoperator_visitfreq(E, ϵ, v)
```

There are several different TE estimators.


##  Transfer entropy from counting-based transfer operator (fast approach)
This method is fast, and can handle a large number of data points. For this example, we'll generate synthetic time series, ``x`` and ``y``, each 20000 data points long.  Then we'll compute the transfer entropy using the visitation
frequency based TE estimator.

```@example s
# Sample a 20000-point orbit from of the logistic3 system, which is
# one of the example systems in CausalityTools.jl.
s = CausalityTools.Systems.logistic3()
orbit = trajectory(s, 20000)
x1, x2 = orbit[:, 1], orbit[:, 2]

# Compute transfer entropy between the 1st and 2nd components.
x1, x2 = orbit[:, 1], orbit[:, 2];
```

### From an embedding with an associated partition at scale `ϵ`

If you are not interested in the internal workings of the transfer entropy
estimator, you can still get a transfer entropy estimate. You still have to
decide on which embedding to use, and how you want to partition the reconstructed state space.

The minimal embedding dimension for the coupled logistic map system is 4 (*exercise: verify this yourself by estimating the box counting dimension of the system using the `DynamicalSystems.jl` package*). Therefore, we'll use a
four-dimensional embedding of the form ``E = {(x2(t+1), x2(t), x1(t), x1(t-1))}``.

Now, check what we need to provide the `transferentropy_visitfreq` estimator.

```
help?> transferentropy_visitfreq
search: transferentropy_visitfreq

  transferentropy_visitfreq(
      E::AbstractEmbedding,
      ϵ::Union{Int, Float64, Vector{Float64}},
      v::TransferEntropy.TEVars) -> Float64

  Using the traditional method of estimation probability distribution by
  visitation frequencies [1], calculate transfer entropy from the embedding E,
  given a discretization scheme controlled by ϵ and information
  v::TEVars about which columns of the embedding to consider for each of the
  marginal distributions. From these marginal distributions, we calculate
  marginal entropies and insert these into the transfer entropy expression.
```

We already know what the embedding will be, so we just need to check what goes
into the `v`, which instructs the estimator which variables of the embedding
goes into which marginal entropy expressions.

```
help?> TEVars
search: TEVars TypeVar

  TEVars(target_future::Vector{Int}
      target_presentpast::Vector{Int}
      source_presentpast::Vector{Int}
      conditioned_presentpast::Vector{Int})

  Which axes of the state space correspond to the future of the target, the
  present/past of the target, the present/past of the source, and the
  present/past of any conditioned variables? Indices correspond to column
  indices of the embedding.points Dataset.

  This information is used by the transfer entropy estimators to ensure the
  marginal distributions are computed correctly.
```

We want to measure TE(x1 -> x2), so we choose `target_future = [1]`, `target_presentpast = [2]` and `source_presentpast = [3, 4]`. We're not conditioning on additional variables, so `conditioned_presentpast = Int[]` is just an empty integer array.

```@example s
# Generate data
s = CausalityTools.Systems.logistic3()
orbit = trajectory(s, 20000)
x1, x2 = orbit[:, 1], orbit[:, 2]

# Which of the time series do we embed?
which_ts = [x1, x2] # these are now indexed 1, 2

# Which variables of the embedding will each of those
# time series correspond to? And what are the lags?
which_pos = [2, 2, 1, 1]
which_lags = [1, 0, 0, -1]
E = embed(which_ts, which_pos, which_lags)

# Construct the TEVars instance.
v = TEVars([1], [2], [3, 4], Int[])

# Subdivide each axis into 15 intervals of equal size, resulting in
# rectangular boxes.
ϵ = 15

# Estimate transfer entropy at scale `ϵ`.
transferentropy_transferoperator_visitfreq(E, ϵ, v)
```
