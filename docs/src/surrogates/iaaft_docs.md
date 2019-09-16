# Surrogate methods

## Iterated amplitude-adjusted Fourier transform (iAAFT) surrogates

```@docs
iaaft
```

## Examples

### Example 1

```@setup iaaft
using Plots, CausalityTools
```

```@example iaaft
npts = 200
ts = sin.(diff(rand(npts + 1)))*0.5 .+ cos.(LinRange(0, 14*pi, npts))
p1 = plot(ts, label = "ts", lc = :black)
p2 = plot(iaaft(ts), label = "iaaft(ts)", xlabel = "Time step")
plot(p1, p2, layout = (2, 1))
ylabel!("Value");

savefig("surr_iaaft.svg"); nothing #hide
```

![](surr_iaaft.svg)

### Example 2

This gif shows iAAFT surrogate realizations for an cyclostationary AR2 (NSAR2) process (`nsar2`) from [1].

![](https://kahaaga.github.io/TimeseriesSurrogates.jl/latest/examples/iaaft_NSAR.gif)

### Example 3

This gif shows iAAFT surrogate realizations for an AR1 process.

![](https://kahaaga.github.io/TimeseriesSurrogates.jl/latest/examples/iaaft_AR1.gif)

## References

1. Lucio et al., Phys. Rev. E *85*, 056202 (2012), after J. Timmer, Phys. Rev. E *58*, 5153 
    (1998). [https://journals.aps.org/pre/abstract/10.1103/PhysRevE.85.056202](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.85.056202)
