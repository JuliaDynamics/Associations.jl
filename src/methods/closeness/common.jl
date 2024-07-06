const HLMS = Union{HMeasure, LMeasure, MMeasure, SMeasure}

function association(measure::HLMS, x::AbstractVector{T}, y::AbstractVector{T}) where T

    (; K, metric, tree_metric, τx, τy, dx, dy, w) = measure
    jsx = ([1 for i = 1:dx]...,)
    τsx = (collect(0:-τx:-(dx-1)*τx)...,)
    jsy = ([1 for i = 1:dy]...,)
    τsy = (collect(0:-τy:-(dy-1)*τy)...,)
    X = genembed(x, τsx, jsx)
    Y = genembed(y, τsy, jsy)
    lX, lY = length(X), length(Y)

    # TODO: cut the last points of the shortest resulting embedding.
    x̂ = lX > lY ? X[1:lY, :] : X
    ŷ = lY > lX ? Y[1:lX, :] : Y
    return association(measure, x̂, ŷ)
end

function association(measure::HLMS, x::AbstractStateSpaceSet{D}, y::AbstractVector{T}) where {D, T}
    (; K, metric, tree_metric, τx, τy, dx, dy, w) = measure

    Y = embed(y, dy, τy)
    X = x[1:length(Y), :]
    return association(measure, X, Y)
end

function association(measure::HLMS, x::AbstractVector{T}, y::AbstractStateSpaceSet{D}) where {D, T}
    (; K, metric, tree_metric, τx, τy, dx, dy, w) = measure

    X = embed(x, dx, τx)
    Y = y[1:length(X), :]
    return association(measure, X, Y)
end
