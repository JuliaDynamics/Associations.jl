Visualise the convergent cross mapping algorithm for a realisation of the [`henon2`](../../example_systems/logistic3.md) system.
The source code for and a description of the `make_ccm_gif` function can be found [here](ccm_gif.md).

```julia
sys_logistic3 = CausalityTools.Systems.logistic3()
tra = trajectory(sys_logistic3, 1000-1)
x, y = tra[:, 1], tra[:, 2]
ts_lengths = 50:10:990

make_ccm_gif(x, y, ts_lengths)
```

![](logistic3.gif)
