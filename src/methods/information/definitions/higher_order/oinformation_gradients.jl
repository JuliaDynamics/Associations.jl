export oinformation_gradients

"""
    oinformation_gradients(est, x₁, x₂, …, xₙ; order = 1) → g

Compute the gradients of the [`OInformation`](@ref) [Scagliarini2023](@cite): how much
the O-information of the inputs changes when single variables (`order = 1`) or pairs of
variables (`order = 2`) are removed. `est` is any estimator of [`OInformation`](@ref),
e.g. `JointProbabilities(OInformation(), disc)` or `GaussianDecomposition(OInformation())`.

## Returns

- `order = 1`: a length-`n` vector with
    ```math
    \\partial_i \\Omega = \\Omega(X) - \\Omega(X_{-i}).
    ```
    ``\\partial_i \\Omega < 0`` means that variable ``i`` takes part in synergistic
    interactions; ``\\partial_i \\Omega > 0`` means that it takes part in redundant ones.
- `order = 2`: a symmetric `n × n` matrix with off-diagonal entries
    ```math
    \\partial^2_{ij} \\Omega = \\Omega(X) - \\Omega(X_{-i}) - \\Omega(X_{-j}) + \\Omega(X_{-ij}),
    ```
    the part of the O-information that needs both ``i`` and ``j``. The diagonal is zero.

Here ``X_{-i}`` denotes all variables except ``X_i``. The O-information of fewer than
three variables is zero. All subset entropies are estimated once and reused.
"""
function oinformation_gradients(est, x...; order::Int = 1)
    if !(est.definition isa OInformation)
        throw(ArgumentError(
            "`oinformation_gradients` needs an estimator of `OInformation`. " *
            "Got an estimator of $(typeof(est.definition))."
        ))
    end
    if order ∉ (1, 2)
        throw(ArgumentError("`order` must be 1 or 2. Got order = $order."))
    end
    n = length(x)
    verify_number_of_inputs_vars(est.definition, n)

    H = subset_entropies(est, x...)
    V = collect(1:n)
    Ω(S) = entropy_combination(est.definition, H, S)
    Ωall = Ω(V)

    if order == 1
        return [Ωall - Ω(setdiff(V, i)) for i in V]
    end

    Ωminus = [Ω(setdiff(V, i)) for i in V]
    g = zeros(n, n)
    for i in 1:n, j in (i + 1):n
        g[i, j] = g[j, i] = Ωall - Ωminus[i] - Ωminus[j] + Ω(setdiff(V, (i, j)))
    end
    return g
end
