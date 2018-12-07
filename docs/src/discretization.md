# Discretization schemes

Across the package, there are many functions that accept the `ϵ` argument.
This is an indication that the underlying algorithm in some way involves a
discretization of a set of points or a delay embedding.  

## Controlling the partitioning
Currently, there are four different ways of partitioning an embedding.
The discretization scheme is controlled by `ϵ`, and the following `ϵ` will work:

* `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
* `ϵ::Float` divides each axis into intervals of size `ϵ`.
* `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
* `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.

## Where are partitions used?

One example of a binning based method is the `transferoperator_grid` estimator
for approximating the transfer (Perron-Frobenius) operator. The
`transferentropy_transferoperator_grid` and `transferentropy_visitfreq` transfer
entropy estimators also both derive probability distributions over
partitioned delay embeddings.
