struct KrishnamurthyL2{B} <: L2DivergenceEstimator
end

function θ̂(k::Function, X; bandwidth = 1)
    N = length(X)
    D = dimension(X)
    f = (1 / (n * (n - 1)))
    g = 1 / (bandwidth ^ D)

    s = 0.0
    for i = 1:N
        xᵢ = X[i]
        s += sum(k((xᵢ, xⱼ) / bandwidth) for (j, xⱼ) in enumerate(X) if j != i) * g
    end

    return f * s
end

function θ̂pq(k::Function, X, Y; bandwidth = 1)
    N = length(X)
    D = dimension(X)
    f = (1 / (n * (n - 1)))
    g = 1 / (bandwidth ^ D)

    s = 0.0
    for i = 1:N
        xᵢ = X[i]
        s += sum(k((xᵢ, xⱼ) / bandwidth) for (j, xⱼ) in enumerate(X) if j != i) * g
    end

    return (1 / (N^2)) * s
end

using StatsBase: sample
function l2_divergence(est::KrishnamurthyL2,
        k::Function,
        X::AbstractDataset{D, T},
        Y::AbstractDataset{D, T}) where {D, T}
    @assert length(X) == length(Y)
    N = length(X)
    firsthalf = sample(1:N, N ÷ 2, replace = false)
    secondhalf = setdiff(1:N, firsthalf)
    θp = @views θ̂(k, X[firsthalf]; bandwidth = est.bandwidth)
    θq = @views θ̂(k, Y[firsthalf]; bandwidth = est.bandwidth)
    θpq = @views θ̂pq(k, X[secondhalf],Y[secondhalf]; bandwidth = est.bandwidth)

    return θp + θq - 2*θpq
end
