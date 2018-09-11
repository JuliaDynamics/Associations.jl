# Quickstart: estimating the transfer operator
## From triangulations

```@example
using CausalityTools # hide
using Plots # hide

ts = [diff(rand(30)) for i = 1:3]
E = invariantize(embed(ts))
triang = triangulate(E)
TO = transferoperator(triang)

# Plot the transfer operator and the invariant distribution
maxprob = min(1, maximum(TO.TO)*1.1) # for plotting
heatmap(TO.TO, clims=(0, maxprob))
xlabel!("Target simplex # (j)")
ylabel!("Source simplex # (i)")
savefig("quick_transferoperator_triang.svg") # hide
nothing # hide
```

![](quick_transferoperator_triang.svg)
