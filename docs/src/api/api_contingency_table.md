
# [Contingency table API](@id contingency_table_api)

To estimate discrete information theoretic quantities that are functions of more than
one variable, we must estimate empirical joint probability mass functions (pmf).
The function [`contingency_matrix`](@ref) accepts an arbitrary number of equal-length
input data and returns the corresponding multidimensional contingency table as a
[`ContingencyMatrix`](@ref). From this table, we can extract the necessary joint and
marginal pmfs for computing any discrete function of multivariate discrete probability
distributions. This is essentially the multivariate analogue of
[`Probabilities`](@ref).

But why would I use a [`ContingencyMatrix`](@ref) instead of some other indirect estimation
method, you may ask. The answer is that [`ContingencyMatrix`](@ref) allows you to
compute *any* of the information theoretic quantities offered in this package for *any*
type of input data. You input data can literally be any hashable type, for example `String`,
`Tuple{Int, String, Int}`, or `YourCustomHashableDataType`.

In the case of numeric data, using a [`ContingencyMatrix`](@ref) is typically a
bit slower than other dedicated estimation procedures.
For example, quantities like discrete Shannon-type [`condmutualinfo`](@ref) are faster to
estimate using a formulation based on sums of four entropies (the H4-principle). This
is faster because we can both utilize the blazingly fast [`StateSpaceSet`](@ref) structure directly,
and we can avoid *explicitly* estimating the entire joint pmf, which demands many
extra calculation steps. Whatever you use in practice depends on your use case and
available estimation methods, but you can always fall back to contingency matrices
for any discrete measure.

```@docs
ContingencyMatrix
contingency_matrix
```

## Utilities

```@docs
marginal_encodings
```
