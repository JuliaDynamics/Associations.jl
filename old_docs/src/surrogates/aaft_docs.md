# Surrogate methods

## Amplitude-adjusted Fourier transform (AAFT) surrogates

```@docs
aaft
```

### Example

```@setup aaft
using Plots, CausalityTools
```

```@example aaft
npts = 200
ts = sin.(diff(rand(npts + 1)))*0.5 .+ cos.(LinRange(0, 14*pi, npts))
p1 = plot(ts, label = "ts", lc = :black)
p2 = plot(aaft(ts), label = "aaft(ts)", xlabel = "Time step")
plot(p1, p2, layout = (2, 1))
ylabel!("Value");

savefig("surr_aaft.svg"); nothing #hide
```

![](surr_aaft.svg)
