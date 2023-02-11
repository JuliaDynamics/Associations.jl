const HLMS = Union{HMeasure, LMeasure, MMeasure, SMeasure}

function estimate(measure::HLMS, x::AbstractVector{T}, y::AbstractVector{T}) where T

    (; K, metric, tree_metric, τx, τy, dx, dy, w) = measure

    X = embed(x, dx, τx)
    Y = embed(y, dy, τy)
    lX, lY = length(X), length(Y)

    # TODO: cut the last points of the shortest resulting embedding.
    x̂ = lX > lY ? X[1:lY, :] : X
    ŷ = lY > lX ? Y[1:lX, :] : Y
    return estimate(measure, x̂, ŷ)
end

function estimate(measure::HLMS, x::AbstractDataset{D}, y::AbstractVector{T}) where {D, T}
    (; K, metric, tree_metric, τx, τy, dx, dy, w) = measure

    Y = embed(y, dy, τy)
    X = x[1:length(Y), :]
    return estimate(measure, X, Y)
end

function estimate(measure::HLMS, x::AbstractVector{T}, y::AbstractDataset{D}) where {D, T}
    (; K, metric, tree_metric, τx, τy, dx, dy, w) = measure

    X = embed(x, dx, τx)
    Y = y[1:length(X), :]
    return estimate(measure, X, Y)
end
