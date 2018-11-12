# Quick start

```@setup embed1
using CausalityTools
using DynamicalSystems
using Plots
```

## Fully customizable embedding

Create a six-dimension embedding lag, with lags `which_lags`, assigning embedding variables by `which_pos`.

```@example embed1
# Initialize a system of three coupled logistic maps and
# create an orbit consisting of 200 points.
s = CausalityTools.Systems.logistic3()
orbit = trajectory(s, 200)

# Use embedding E = {(x(t+3), x(t+2), x(t), y(t), y(t-2), z(t))}
which_pos = [1, 1, 1, 2, 2, 3]
which_lags = [3, 2, 0, 0, -2, 0]
E = embed(orbit, which_pos, which_lags)
scatter3d(E.points[1, :], E.points[4, :], E.points[6, :], ms = 1.5,
        xlabel = "rn1(t+2)", ylabel = "rn2(t)", zlabel = "rn2(t-6)",
        legend = false) # hide
savefig("embed_quickex2.svg") # hide
nothing # hide
```

![](embed_quickex2.svg)


## Wrap pre-computed embedding

Here's how you would wrap a pre-computed embedding. 

```@example embed1
# Initialize a system of three coupled logistic maps and
# create an orbit consisting of 200 points.
s = CausalityTools.Systems.logistic3()
orbit = trajectory(s, 200)

# Embed the unmodified dataset, assuming we want the embedding
# E = {(x(t), y(t), z(t))}
E = embed(orbit)
scatter3d(E.points[1, :], E.points[2, :], E.points[3, :], ms = 1.5,
        xlabel = "x(t)", ylabel = "y(t)", zlabel = "z(t)",
        legend = false) # hide
savefig("embed_quickex1.svg") # hide
nothing # hide
```

![](embed_quickex1.svg)
