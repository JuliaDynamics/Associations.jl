import TransferEntropy: te_embed

"""
    NearestNeighbourMITest(; k::Int = 1, l::Int = 1, m::Int = 1, n::Int = 1; τ::Int = 1,
        estimator::NearestNeighbourMI = NearestNeighbourMI(b = 2, k1 = 2, k2 = 3, metric = Chebyshev()),
        ηs::Union{Int, AbstractVector{Int}})

The parameters for a transfer entropy test using the `NearestNeighbourMI` estimator that estimates 
mutual information terms using the counting of nearest neighbours.

## Mandatory keyword arguments

- **`ηs`**: The prediction lags (that gos into the ``T_{f}`` component of the embedding).

## Optional keyword arguments

- **`k::Int`**: The dimension of the ``T_{f}`` component of the embedding. 

- **`l::Int`**: The dimension of the ``T_{pp}`` component of the embedding. 

- **`m::Int`**: The dimension of the ``S_{pp}`` component of the embedding. 

- **`n::Int`**: The dimension of the ``C_{pp}`` component of the embedding. 

- **`τ::Int`**: The embedding lag. Default is `τ = 1`.

- **`metric`**: The distance metric to use for mutual information estimation.

- **`k1`**: The number of nearest neighbours for the highest-dimensional mutual
    information estimate. To minimize bias, choose ``k_1 < k_2`` if
    ``min(k_1, k_2) < 10`` (see fig. 16 in [1]). Beyond dimension 5, choosing
    ``k_1 = k_2`` results in fairly low bias, and a low number of nearest
    neighbours, say `k1 = k2 = 4`, will suffice.

- **`k2`**: The number of nearest neighbours for the lowest-dimensional mutual
    information estimate. To minimize bias, choose ``k_1 < k_2`` if
    if ``min(k_1, k_2) < 10`` (see fig. 16 in [1]). Beyond dimension 5, choosing
    ``k_1 = k_2`` results in fairly low bias, and a low number of nearest
    neighbours, say `k1 = k2 = 4`, will suffice.
"""
Base.@kwdef mutable struct NearestNeighbourMITest{N} <: TransferEntropyCausalityTest{N}
    
    """ The delay reconstruction parameter k (controls dimension of ``T_{f}`` component of embedding). """
    k::Int = 1

    """ The delay reconstruction parameter l (controls dimension of ``T_{pp}`` component of embedding). """
    l::Int = 1

    """ The delay reconstruction parameter m (controls dimension of ``S_{pp}`` component of embedding). """
    m::Int = 1

    """ The delay reconstruction parameter n (controls dimension of ``C_{pp}`` component of embedding). """
    n::Int = 1

    """ The delay reconstruction lag for the ``T_{pp}`` component of the embedding. """
    τ::Int = 1

    """ The transfer entropy estimator. """
    estimator::NearestNeighbourMI = NearestNeighbourMI(b = 2, k1 = 2, k2 = 3, metric = Chebyshev())

    """ The prediction lags"""
    ηs::Union{Int, AbstractVector{Int}}

    function NearestNeighbourMITest(k::Int, l::Int, m::Int, n::Int, τ::Int, 
            estimator::NearestNeighbourMI, ηs)
            
        N = length(ηs) # length of return vector when used with `causality`
        return new{N}(k, l, m, n, τ, estimator, ηs)
    end
end

function causality(source::AbstractVector{T}, target::AbstractVector{T}, test::NearestNeighbourMITest{N}) where {N, T}
    
    length(source) == length(target) || error("Input vectors must be of same length")
    
    τ = test.τ
    k = test.k
    l = test.l
    m = test.l
    
    tes = zeros(T, N)
    
    for (i, η) in enumerate(test.ηs)
        pts, vars = te_embed(source, target, k, l, m; η = η, τ = τ);
        tes[i] = transferentropy(pts, vars, test.estimator)
    end
    
    return tes
end

export NearestNeighbourMITest