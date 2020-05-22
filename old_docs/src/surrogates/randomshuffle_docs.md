# Surrogate methods

## Random shuffle surrogates

```@docs
randomshuffle
```

### Example

```@setup randomshuffle
using Plots
using CausalityTools
```

```@example randomshuffle
npts = 100
ts = sin.(diff(rand(npts + 1)))*0.5 .+ cos.(LinRange(0, 8*pi, npts))
p1 = plot(ts, label = "ts", lc = :black)
p2 = plot(randomshuffle(ts), label = "randomshuffle(ts)", xlabel = "Time step")
plot(p1, p2, layout = (2, 1))
ylabel!("Value");

savefig("surr_randomshuffle.svg"); nothing #hide
```

![](surr_randomshuffle.svg)
