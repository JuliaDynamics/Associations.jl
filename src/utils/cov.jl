import Statistics.cov
using Statistics: mean
using StateSpaceSets: AbstractStateSpaceSet
using StaticArrays: @MMatrix, @MVector, SMatrix, SVector

export fastcov

# Non-allocating and more than twice as fast as writing a wrapper
# `f(x) = Statistics.cov(Matrix(x))`.
# Also accepts SubStateSpaceSets, so we can use views on neighbor points.
# These functions return StaticArrays.
fastcov(x::AbstractStateSpaceSet) = fastcov(x.data)
fastmean_and_cov(x::AbstractStateSpaceSet) = fastmean_and_cov(x.data)

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

# `fastcor(x)` is twice as fast as `cor(Matrix(x)` and non-allocating.
export fastcor
fastcor(x::AbstractStateSpaceSet) = fastcor(x.data)
function fastcor(x::Vector{SVector{D, T}}) where {D, T}
    N = length(x)
    μ, Σ = fastmean_and_cov(x)
    σ = std(x)
    C = @MMatrix zeros(D, D)
    for j in 1:D
        for i in 1:D
            C[i, j] = Σ[i, j] / (σ[i] * σ[j])
        end
    end
    return SMatrix{D, D}(C)
end
