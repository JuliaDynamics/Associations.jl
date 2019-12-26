import TransferEntropy: te_embed

"""
    NearestNeighbourMITest(; k::Int = 1, l::Int = 1, m::Int = 1, n::Int = 1; τ::Int = 1,
        estimator::NearestNeighbourMI = NearestNeighbourMI(b = 2, k1 = 2, k2 = 3, metric = Chebyshev()),
        ηs::Union{Int, AbstractVector{Int}})

A transfer entropy test that estimates mutual information terms by counting nearest neighbours. 
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