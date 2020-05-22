# Delay embeddings (state space reconstructions)

The concept of state space reconstructions (SSR), or delay embeddings, is crucial for
working with dynamical systems. Delay embeddings become important when analysing the 
dynamics of real world data, where the underlying, governing equations of the system are 
not known, because they can be used to reconstruct the dynamics of the underlying system 
directly from observations of only a subset of the entire dynamical system.

Delay embeddings are required for many of the algorithms in `CausalityTools`. For basic
usage, you may not need to think about embeddings at all. It is highly recommended, though,
to try and understand the algorithms. For that, you need to manually perform delay embeddings.

## Customized embeddings
To construct embeddings from data, you can use the `customembed(data, positions, lags)`
function.

The first argument `data` is the dataset, and may a matrix where either the columns or
rows are dynamical variables (the largest dimension will be interpreted as the 
(time, depth, position, etc) index), a vector of dynamical variables, or a `Dataset` instance
from [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/latest/).

The second argument is a vector of positions, indicating which dynamical variables are
assigned to which coordinate axis in the embedding. The length of this vector dictates
the final dimension of the reconstructed space. For example, `[1, 1, 2, 5]` means that
the first dynamical variable should be assigned to both the first and second coordinate axis
of the embedding, the second dynamical variables goes in the third coordinate axis, and
that fifth dynamical variable is assigned to the fourth coordinate axis.

The third argument is a vector of embedding lags, one for each coordinate axis of the
embedding. If `[1, 1, 2, 5]` are the positions, then embeddings lags `[0, 2, 0, -3]`
indicate that embedding axis ``\#``1 (dynamical variable ``\#``1) should not be lagged, that
embedding axis ``\#``2 (also constructed from dynamical variable ``\#``1) should have a positive
lag of `2`, that axis ``\#``3 (dynamical variable ``\#``2) should be non-lagged, and that
axis ``\#``4 (dynamical variables ``\#``5) should have a negative lag of `-3`.

Here are some more examples. Assume that `data` holds the observations of the dynamical
variables ``x_1, x_2, \ldots, x_n``.


| Embedding | Code |
|-----------|------|
|``\left{ (x_1(t+3), x_1(t+1), x_1(t), x_2(t), x_2(t-3)) \right}`` | `customembed(data, [1, 1, 1, 2, 2], [3, 1, 0, 0, -3]` |
| ``\left{ (x_1(t), x_2(t), x_3(t)) \right}`` | `customembed(data, [1, 2, 3], [0, 0, 0]` |
| ``\left{ (x_4(t+5), x_3(t-2), x_6(t-8), x_1(t+1)) \right}`` | `customembed(data, [4, 3, 6, 1], [5, -2, -8, 1])` |


This way of constructing embeddings is very useful when computing, for example, transfer
entropy, because it allows the user to have full control over which embedding axes
goes into which marginal entropy expressions (thus, reinforcing
understanding the algorithm, rather than blindly accepting that it works).


Embeddings may also be constructed using the `embed` or `reconstruct` functions from
[DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/latest/).
However, for many of the algorithms we need a bit more flexibility when constructing our
state space reconstructions, which is what `customembed` is for.
