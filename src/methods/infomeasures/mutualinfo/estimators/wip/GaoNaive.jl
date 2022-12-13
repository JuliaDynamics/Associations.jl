# export GaoNaiveMI
# using Neighborhood: bulksearch
# using Neighborhood: Euclidean, KDTree, NeighborNumber, Theiler

# Base.@kwdef struct GaoNaiveMI <: MutualInformationEstimator
#     k::Int = 1
#     w::Int = 0
# end

# function mutualinfo(e::Renyi, est::GaoNaiveMI, x::Vector_or_Dataset...)
#     e.q == 1 || throw(ArgumentError(
#         "Renyi entropy with q = $(e.q) not implemented for $(typeof(est)) estimators"
#     ))
#     (; k, w, base) = est

#     joint = Dataset(x...)
#     D = dimension(joint)
#     N = length(joint)
#     M = length(x)

#     f_joint = (k  * gamma(D / 2 + 1)) / ( (N - 1) * π^(D / 2))
#     tree = KDTree(joint, Euclidean())
#     ds_joint = last.(bulksearch(tree, joint, NeighborNumber(k), Theiler(w))[2])

#     ds = Vector{Vector{Float64}}(undef, M)
#     fs = Vector{Float64}(undef, M)
#     for (i, xᵢ) in enumerate(x)
#         Xᵢ = Dataset(xᵢ)
#         Dᵢ = dimension(xᵢ)
#         fᵢ = (k  * gamma(Dᵢ / 2 + 1)) / ( (N - 1) * π^(D / 2))
#         ds[i] = last.(bulksearch(tree, Xᵢ, NeighborNumber(k), Theiler(w))[2])
#     end

#     mi = 0.0

#     for i in 1:N
#         mi += log(ds_joint[i]) / prod(ds[j][i] for j = 1:M)
#     end

#     mi /= N

#     return mi / log(base, ℯ)

# end
