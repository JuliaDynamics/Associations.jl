# Invariant measures over a partitioned state space


```@setup s
using CausalityTools
using Plots
using DynamicalSystems
```

Obtaining an invariant measure over a state space at a given resolution is simple.

A partition at a given resolution ``ϵ``  subdivides the embedding into ``n`` states. Compute the transfer matrix ``Tϵ``, giving the transition probabilities between the states ``s_i`` of the partition. Choose an initial distribution ``p_0(x_i)``, represented as a row vector, where ``i = 1, 2, \ldots, n``. Obtaining the invariant measure the partition elements is then obtained by repeated application of the transfer operator to this distribution. In other words, ``p_{inv}(xᵢ) = p_0(x_i)Tϵ^{n}``. Numerically, we stop the iteration when
the relative change between consecutive distributions is below some threshold.

## Computing the invariant measure
To compute the invariant measure, we follow the exact same procedure as when generating the transfer matrix. However, we add the additional step of finding the left eigenvector of ``Tϵ``. This can be done using the function `left_eigenvector(TO::AbstractTransferOperator)`.

## Short time series (triangulation approach)
```@example s
ts = [diff(rand(30)) for i = 1:3]
E = invariantize(embed(ts))
triang = triangulate(E)
TO = transferoperator(triang)
invdist = left_eigenvector(TO).dist

# Plot the transfer operator and the invariant distribution
maxprob = min(1, maximum(TO.TO)*1.1) # for plotting
p1 = heatmap(TO.TO, clims=(0, maxprob))
xlabel!(p1, "Target simplex # (j)")
ylabel!(p1, "Source simplex # (i)")
p2 = bar(invdist, legend = false)
xlabel!(p2, "Source simplex #")
ylabel!(p2, "Invariant measure")
l = @layout [a{1.0w}; b{1.0w}]
plot(p1, p2, layout = l)
```

The heapmap visualizes the transfer matrix, which gives the transition probabilities between all pairs of states ``(s_i, s_j)``. The bar plot shows the invariant measure (i.e. the long-term visitation frequencies) over the states ``s_i``.


## Long time series (direct approach)
```@example s
using CausalityTools, Plots
ts = [diff(rand(5000)) for i = 1:4]
E = invariantize(embed(ts))
ϵ = 20
TO = transferoperator(E, ϵ)
invdist = left_eigenvector(TO).dist

# Plot the transfer operator and the invariant distribution
maxprob = min(1, maximum(TO.TO)*1.1) # for plotting
p1 = heatmap(TO.TO, clims=(0, maxprob))
xlabel!(p1, "Source partition box # (j)")
ylabel!(p1, "Target partition box # (i)")

p2 = bar(invdist, legend = false)
xlabel!(p2, "Source partition box #")
ylabel!(p2, "Invariant measure")
plot(p1, p2, layout = @layout [a{1.0w}; b{1.0w}])
```

Again, the heatmap shows the transition probabilities between all states ``(s_i, s_j)``, while
the bar plot shows the invariant distribution over the states ``s_i`` of the partition.
