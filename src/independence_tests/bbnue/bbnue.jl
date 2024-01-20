# TODO: test how this works with 
export BootstrapBasedNonUniformEmbeddingTest
export MostInformativeForwardSearch
export forward_search

"""
    BootstrapBasedNonUniformEmbeddingTest(;
        measure_pairwise = MIShannon();
        measure_cond
    )

"""
Base.@kwdef struct BootstrapBasedNonUniformEmbeddingTest{M} <: IndependenceTest{M}
    measure_pairwise::M = MIShannon()
    measure_cond::M = CMIShannon()
    max_dim::Int = 2
    max_τ::Int = 3
    include_zerolag::Bool = false
end

function prepare_embedding_lags_and_indices(test, x, y, 𝒵::AbstractStateSpaceSet)
    min_τ = test.include_zerolag ? 0 : 1
    τs_x = min_τ:test.max_τ
    τs_y = min_τ:test.max_τ
    τs = min_τ:test.max_τ
    n_lags = length(τs )
    τs_z = vcat([collect(τs) for i = 3:3+dimension(𝒵)-1]...,)

    # Variables are enumerated according to their input order.
    # e.g. x = 1, y = 2, z[:, 1] = 3, z[:, 2] = 4, and so on.
    js_x = repeat([1], n_lags)
    js_y = repeat([2], n_lags)
    js_z = vcat([repeat([i], n_lags) for i = 3:3+dimension(𝒵)-1]...,)

    τs = [τs_x; τs_y; τs_z]
    js = [js_x; js_y; js_z]
    return τs, js
end

function joint_embedding(test::BootstrapBasedNonUniformEmbeddingTest, x, y, 𝒵)
    joint_dataset = [x y Matrix(z)]
    τs, js = prepare_embedding_lags_and_indices(test, x, y, 𝒵)

    # Construct the candidate set.
    𝒞 = genembed(Dataset(joint_dataset), τs, js)
end

function nue(test::BootstrapBasedNonUniformEmbeddingTest, 
        x::AbstractVector{T}, 
        y::AbstractVector{T}, 
        𝒵::AbstractDataset{D, T}) where {D, T}

    𝒮 = Vector{T}[]
    
end

# Notation:
# 𝒮: selected variables
# 𝒞: candidate variables

function k_forward_search!(test::BootstrapBasedNonUniformEmbeddingTest, 𝒞, 𝒮, k::Int)
    if k == 0
    end
end

function forward_search!(test::BootstrapBasedNonUniformEmbeddingTest, 𝒞, 𝒮)

end

# TODO: provide estimator, not measure
"""
    MostInformativeForwardSearch(;
        measure_pairwise = MIShannon(),
        measure_cond = CMIShannon(),
        n_reps = 100,
        α = 0.05
    )

The forward search part of the BBNUE algorithm from Baboukani et al. (2020).
"""
Base.@kwdef struct MostInformativeForwardSearch{MP, MC}
    measure_pairwise::MP = GaussianMI(MIShannon())
    measure_cond::MC = GaussianCMI(CMIShannon())
    n_reps = 100
    α = 0.05
end

abstract type StoppingCriterion end
struct BootstrapBasedStopping <: StoppingCriterion end 

"""
    forward_search!(test::MostInformativeSearch, target, 𝒞, 𝒮)

Given a `target` variable, a set of candidate variables `𝒞` and a set of 
(initially empty) selected variables, apply the `test` and do the search.

`idxs_𝒞` and `idxs_𝒮` keep track of the indices of the selected/candidate variables,
in terms of their input order.
"""
function forward_search(search::MostInformativeForwardSearch, target, candidates)
    # Initialize empty set of selected variables
    𝒮 = eltype(candidates)[]

    # Make a copy of the input candidates, so that we don't remove any of the original
    # data.
    𝒞 = StateSpaceSet.(deepcopy(candidates))

    k = 0
    terminate = false
    idxs_𝒞 = collect(1:length(𝒞))
    idxs_𝒮 = Int[]
    @show idxs_𝒞, idxs_𝒮

    while !terminate
        kth_forward_search!(search, target, 𝒞, 𝒮, idxs_𝒞, idxs_𝒮, k)
        # Shuffle target and shuffle the *last* selected variable
        k += 1 
        if k >= length(candidates)
            terminate = true
        end
    end

    return idxs_𝒞, idxs_𝒮
end

function kth_forward_search!(search::MostInformativeForwardSearch, target, 𝒞, 𝒮, idxs_𝒞, idxs_𝒮, k)
    if k == 0
        forward_search_pairwise!(search, target, 𝒞, 𝒮, idxs_𝒞, idxs_𝒮)
    else
        forward_search_conditional!(search, target, 𝒞, 𝒮, idxs_𝒞, idxs_𝒮)  
    end
end

# Perform a single forward search with the pairwise measure
function forward_search_pairwise!(search::MostInformativeForwardSearch, target, 
        𝒞, 𝒮, idxs_𝒞, idxs_𝒮)
    ℳ = search.measure_pairwise
    associations = [CausalityTools.estimate(ℳ, target, 𝒞[i]) for i in eachindex(𝒞)]
    @show associations

    𝒜, idx_most_informative_variable = findmax(associations)

    update_variables!(idx_most_informative_variable, 𝒞, 𝒮, idxs_𝒞, idxs_𝒮)
end

# Perform a single forward search with the conditional measure
function forward_search_conditional!(search::MostInformativeForwardSearch, target, 
        𝒞, 𝒮, idxs_𝒞, idxs_𝒮)
    ℳ = search.measure_cond
    s_concatenated = StateSpaceSet(𝒮...)

    associations = [CausalityTools.estimate(ℳ, target, 𝒞[i], s_concatenated) for i in eachindex(𝒞)]
    @show associations
    𝒜, idx_most_informative_variable = findmax(associations)

    return update_variables!(idx_most_informative_variable, 𝒞, 𝒮, idxs_𝒞, idxs_𝒮)
end

function update_variables!(idx_most_informative_variable, 𝒞, 𝒮, idxs_𝒞, idxs_𝒮)
    push!(𝒮, 𝒞[idx_most_informative_variable])
    push!(idxs_𝒮, idx_most_informative_variable)
    deleteat!(𝒞, idx_most_informative_variable)
    deleteat!(idxs_𝒞, idx_most_informative_variable)
    @show idxs_𝒞, idxs_𝒮
end