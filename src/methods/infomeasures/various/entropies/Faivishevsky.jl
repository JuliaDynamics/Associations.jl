export Faivishevsky
#https://proceedings.neurips.cc/paper/2008/file/3dc4876f3f08201c7c76cb71fa1da439-Paper.pdf
Base.@kwdef struct Faivishevsky{M} <: DifferentialEntropyEstimator
    k::Int = 1 # todo: remove. it isn't used.
    w::Int = 0
    metric::M = Euclidean()
end
import ComplexityMeasures.ball_volume
using Neighborhood: search
function entropy(e::Shannon, est::Faivishevsky, x::AbstractStateSpaceSet{D}) where D
    (; k, w, metric) = est
    N = length(x)
    tree = KDTree(x, metric)
    idxs, ϵs  = bulksearch(tree, x, NeighborNumber(N-1), Theiler(w))

    f = 0.0
    for k = 1:N-1 # loop over neighbor numbers
        f -= digamma(k)
        c = 0.0
        for i = 1:N # loop over points
            c += D/N * log(ϵs[i][k])
        end
    end
    f *= 1 / (N - 1)
    # The unit is nats
    h = digamma(N) + ball_volume(D) + f
    return _convert_logunit(h, ℯ, e.base)
end
