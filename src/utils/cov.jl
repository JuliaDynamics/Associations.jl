import Statistics.cov
using Statistics: mean
using StateSpaceSets: AbstractDataset
using StaticArrays: @MMatrix, @MVector, SMatrix

export fastcov

# Non-allocating and more than twice as fast as writing a wrapper
# `f(x) = Statistics.cov(Matrix(x))`.
# Also accepts SubDatasets, so we can use views on neighbor points.
# These functions return StaticArrays.
fastcov(x::AbstractDataset) = fastcov(x.data)
fastmean_and_cov(x::AbstractDataset) = fastmean_and_cov(x.data)

function fastcov(x̄, x::Vector{SVector{D, T}}) where {D, T}
    T <: AbstractFloat || error("Need `eltype(x[i]) <: AbstractFloat` ∀ i ∈ 1:length(x). Got `eltype(x[i])=$(eltype(first(x)))`")
    N = length(x) - 1
    C = @MMatrix zeros(D, D)
    x̄ = mean(x)
    Δx = @MVector zeros(D)
    @inbounds for xᵢ in x
        Δx .= xᵢ - x̄
        C .+= Δx * transpose(Δx)
    end
    C ./= N
    return SMatrix{D, D}(C)
end
# So we don't have to compute the mean twice at every iteration.
function fastcov(x::Vector{SVector{D, T}}) where {D, T}
    T <: AbstractFloat || error("Need `eltype(x[i]) <: AbstractFloat` ∀ i ∈ 1:length(x). Got `eltype(x[i])=$(eltype(first(x)))`")

    μ = mean(x)
    fastcov(μ, x)
end
function fastmean_and_cov(x::Vector{SVector{D, T}}) where {D, T}
    μ = mean(x)
    Σ = fastcov(μ, x)
    return μ, Σ
end
