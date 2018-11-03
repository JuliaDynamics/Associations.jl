# Reducing bias
Many information theoretic calculations requires a coarse graining of the
state/phase space. A priori, there is no way of knowing what the "correct"
partition is. To reduce bias, it is recommended to compute transfer entropy
over a set of different partitions ranging in coarseness, then consider an
average transfer entropy across partitions of varying bin sizes.

```@setup s
using CausalityTools
using Plots
```
