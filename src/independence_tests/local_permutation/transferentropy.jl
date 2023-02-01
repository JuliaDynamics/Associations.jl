
function independence(test::LocalPermutationTest{<:TransferEntropy}, x::AbstractVector...)
    throw(ArgumentError("`independence` not yet implemented for test `$(typeof(test))`"))
end
# function independence(test::LocalPermutationTest{<:TransferEntropy{<:E, <:EmbeddingTypes}}, x::AbstractVector...) where E
#     (; measure, est, rng, kperm, nsurr) = test
#     Ŝ, T⁺, S, T = marginals_and_surrogenerator(measure.embedding, surrogate, x...)
#     @assert length(T⁺) == length(S) == length(T)
#     N = length(x)
#     cmi = te_to_cmi(measure)

#     Î = estimate(cmi, est, T⁺, S, T)
#     tree_z = KDTree(Z, Chebyshev())
#     idxs_z = bulkisearch(tree_z, Z, NeighborNumber(kperm), Theiler(0))
#     𝒩 = MVector{kperm, Int16}.(idxs_z) # A statically sized copy
#     n̂ = collect(1:N)
#     X̂ = deepcopy(X)
#     𝒰 = zeros(Int, N) # used indices
#     Îs = zeros(nsurr)
#     for b in 1:nsurr
#         shuffle_neighbor_indices!(𝒩, rng)
#         # By re-filling, we avoid allocating extra vector for each surr. By filling with
#         # zeros, we make sure that the while loop below isn't affected.
#         𝒰 .= 0
#         Π = new_permutation!(n̂, rng)
#         for i in Π # for every point xᵢ.
#             𝒩ᵢ = 𝒩[i] # shuffled neighbors to xᵢ, in terms of z
#             j = first(𝒩ᵢ)
#             m = 1
#             while j ∈ 𝒰 && m < kperm
#                 m += 1
#                 j = 𝒩ᵢ[m]
#             end
#             𝒰[i] = j
#             push!(𝒰, j)
#             X̂.data[i] = X.data[j]
#         end
#         b == 1 && @show X̂
#         @show measure
#         Îs[b] = estimate(cmi, est, X̂, Y, Z)
#     end
#     p = count(Î .<= Îs) / nsurr

#     return LocalPermutationTestTest(Î, Îs, p, nsurr)
# end
