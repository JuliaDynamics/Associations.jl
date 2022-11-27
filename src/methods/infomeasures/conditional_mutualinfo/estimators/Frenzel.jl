
using Neighborhood: bulkisearch, inrangecount
using Neighborhood: Theiler, NeighborNumber, KDTree, Chebyshev
using SpecialFunctions: digamma

export Frenzel

Base.@kwdef struct Frenzel{MJ, MM} <: ConditionalMutualInformationEstimator
    k::Int = 1
    w::Int = 0
    metric_joint::MJ = Chebyshev()
    metric_marginals::MM = Chebyshev()
end


function estimate(infomeasure::CMI{Nothing}, e::Renyi, est::Frenzel, X, Y, Z)
    e.q ≈ 1 || throw(ArgumentError(
        "Renyi entropy with q = $(e.q) not implemented for $(typeof(est)) estimators"
    ))
    (; k, w, metric_joint, metric_marginals) = est
    @assert length(X) == length(Y) == length(Z)
    N = length(X)
    # Ensures that vector-valued inputs are converted to datasets, so that
    # building the marginal/joint spaces and neighbor searches are fast.
    X = Dataset(X)
    Y = Dataset(Y)
    Z = Dataset(Z)
    joint = Dataset(X, Y, Z)
    XZ = Dataset(X, Z)
    YZ = Dataset(Y, Z)

    tree_joint = KDTree(joint, metric_joint)
    ds_joint = last.(bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(w))[2])
    tree_xz = KDTree(XZ, metric_marginals)
    tree_yz = KDTree(YZ, metric_marginals)
    tree_z = KDTree(Z, metric_marginals)

    cmi = time_averaged_harmonic_numbers(tree_xz, tree_yz, tree_z, XZ, YZ, Z, ds_joint, N) -
        harmonic_number(k - 1)

    return cmi / log(e.base, ℯ)
end

function time_averaged_harmonic_numbers(tree_xz, tree_yz, tree_z, XZ, YZ, Z, ds_joint, N)
    mean_dgs = 0.0
    for (i, dᵢ) in enumerate(ds_joint)
        # Usually, we subtract 1 because inrangecount includes the point itself,
        # but we'll have to add it again inside the digamma, so just skip it.
        nxz = inrangecount(tree_xz, XZ[i], dᵢ)
        nxy = inrangecount(tree_yz, YZ[i], dᵢ)
        nz = inrangecount(tree_z, Z[i], dᵢ)
        mean_dgs += harmonic_number(nxz) +
            harmonic_number(nxy) -
            harmonic_number(nz)
    end

    return mean_dgs / N
end

harmonic_number(n::Int) = -sum(1/i for i = 1:n)
# Example from their paper.
# using Distributions: Normal
# using Statistics: cor
# using LinearAlgebra: eigvals
# function sunspot_model(n::Int; ϵ, base = 2,
#         α₁ = 1.90694, α₂ = -0.98751,
#         β₁ = 0.78512, β₂ = -0.40662)
#     Na1 = Normal(0, 1)
#     Na2 = Normal(0, 1)
#     a1 = rand(Na1, n)
#     a2 = rand(Na2, n)
#     z1 = zeros(n)
#     z2 = zeros(n)
#     z1[1] = rand(Na1); z1[2] = rand(Na1)
#     z2[1] = rand(Na1); z2[2] = rand(Na1)
#     for i = 3:n
#         z1[i] = α₁*z1[i-1] + α₂*z1[i-2] + a1[i] - β₁*a1[i-1] - β₂*a1[i-2]
#         z2[i] = α₁*(ϵ*z1[i-1]+(1-ϵ)*z2[i-1]) +
#             α₂*z2[i-2] + a2[i] - β₁*a2[i-1] - β₂*a2[i-2]
#     end
#     Z = [z1[2:end] z2[2:end] circshift(z2, -1)[2:end]]
#     cmi = -0.5*sum(log.(eigvals(cor(Z))))
#     return z1[2:end], z2[2:end], circshift(z2, -1)[2:end], cmi / log(base, ℯ)
# end
