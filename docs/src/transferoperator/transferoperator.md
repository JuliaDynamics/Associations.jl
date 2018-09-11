# Estimating the transfer operator

For short time series, the most reliable estimates of the transfer operator are obtained by using a triangulation of the state space as our partition. This approach is computationally costly because it has to compute N-dimensional simplex intersections. However, it gives robust estimates of ergodic transition probabilities down to as little as a few hundred points.

For longer time series, we use a rectangular grid to discretize the embedding.
This approach is must faster, because no intersections have to be explicitly computed.

```@setup s
using CausalityTools
using Plots
using DynamicalSystems
```

## Short time series
Embed some random short time series (n = 30), make sure that the embedding is invariant, triangulate the embedding, and estimate the transfer operator on the
states (partition elements) of the discretized state space. We then compute the invariant measure over the states from the transfer matrix.

```@example s
ts = [diff(rand(30)) for i = 1:3]
E = invariantize(embed(ts))
triang = triangulate(E)
TO = transferoperator(triang)

# Plot the transfer operator and the invariant distribution
maxprob = min(1, maximum(TO.TO)*1.1) # for plotting
heatmap(TO.TO, clims=(0, maxprob))
xlabel!("Target simplex # (j)")
ylabel!("Source simplex # (i)")
savefig("transferoperator-short-ex.svg"); nothing #hide
```

![](transferoperator-short-ex.svg)

The heapmap visualizes the transition probabilities between the elements of the partition (i.e. the simplices). Each row in the transfer matrix sums to 1, completely accounting for the possible transitions the i-th state can make. The transfer matrix is Markov, so the `(i,j)`th element of it gives the probability of the orbit going from state `i` to state `j` in the next time step.


## Long time series
For the longer time series, we use the `transferoperator(E::AbstractEmbedding, ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}})` estimator, where the partition is constructed according to `ϵ`. It first identify which bins of the partition are visited by which points of the orbit, then uses that information to generate the transfer matrix. Computationally, the procedure is almost identical to what we would do for the triangulation approach.

Embed some random long time series (n = 10000), and make sure that the embedding is invariant. Then decide on how to partition the state space. We set `ϵ = 5`,
which divides each axis of the embedding into 5 intervals of equal length and
estimate the transfer matrix, from which we gain the invariant measure over
the partition elements.

```@example s
ts = [diff(rand(10000)) for i = 1:4]
E = embed(ts)
ϵ = 5
TO = transferoperator(E, ϵ)

# Plot the transfer operator and the invariant distribution
maxprob = min(1, maximum(TO.TO)*1.1) # for plotting
heatmap(TO.TO, clims=(0, maxprob))
xlabel!("Source partition box # (j)")
ylabel!("Target partition box # (i)")
savefig("transferoperator-long-ex.svg"); nothing #hide
```

![](transferoperator-long-ex.svg)

Using `ϵ = 5` as the binning scheme, the orbit visits just over 500 unique states. The heatmap again shows the transition probabilities, according to the transfer operator, between all states
``(s_i, s_j)``.
