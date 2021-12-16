import TransferEntropy: process_input, embed_candidate_variables, optim_te
using DelayEmbeddings, Statistics

# Usaully, we use all lags from startlag:-\tau_max to construct variables. In some situations,
# we may want to exclude som of those variables. 
function filtered_τs(τs::AbstractVector{Int}, js::AbstractVector{Int}, τexclude::Int)
    [τ for τ in τs if abs(τ) != abs.(τexclude)]
end

function filtered_js(τs::AbstractVector{Int}, js::AbstractVector{Int}, τexclude::Int)
    [j for (τ, j) in zip(τs, js) if abs(τ) != abs.(τexclude)]
end

export pa_candidate_variables

function get_delay_τs(source, target; maxlag::Union{Int, Float64} = 0.05)
# Ensure all time series are of the same length.
    Ls = [length.(source); length.(target)]
    @assert all(Ls .== maximum(Ls))

    if maxlag isa Int
        delay_τs = 1:maxlag
    else
        delay_τs = 1:ceil(Int, maximum(Ls)*maxlag)
    end
    return delay_τs
end

"""
    search_τs(source, target, cond; 
        maxlag::Union{Int, Float64} = 0.05, τcutoff::Int = 100) → delay_τs

Decide delay range search among when searching for maximum embedding delays, based on input 
time series length(s).
"""
function search_τs(source, target, cond; 
        maxlag::Union{Int, Float64} = 0.05, τcutoff::Int = 100)

    # Ensure all input time series are of the same length.
    Ls = [length.(source); length.(target); length.(cond)]
    @assert all(Ls .== maximum(Ls))
    
    if maxlag isa Int
        delay_τs = 1:minimum(maxlag, cutoff)
    else
        delay_τs = 1:ceil(Int, maximum(Ls)*maxlag)
    end
    return delay_τs
end

"""
    bbnue_embedding_lags(include_instantaneous, τmax, exclude)

Generate embedding lags for predictive asymmetry analysis using the bbnue algorithm. 
"""
function bbnue_embedding_lags(include_instantaneous, τmax, exclude)
    @assert τmax > 0
    startlag = include_instantaneous ? 0 : -1
    [τ for τ in startlag:-1:-τmax if abs(τ) ∉ abs.(exclude)]
end

"""
    pa_candidate_variables(
        source::Vector{AbstractVector}, 
        target::Vector{AbstractVector}, 
        [cond::Vector{AbstractVector}];
        k::Int = 1, include_instantaneous = true,
        τexclude::Union{Int, Nothing} = nothing,
        maxlag::Union{Int, Float64} = 0.05
    ) → ([τs_source, τs_target, τs_cond, ks_T⁺, ks_T⁻], [js_source, js_target, js_cond, js_T⁺, js_T⁻])

Construct candidate variables from input time series. `source` is a vector of equal-length time series
assumed to represent the putative source process. `target` and `cond` are the same, but contains time series
of the target process and of the conditional processes, respectively. `k` is the desired prediction lag. 

If `include_instantaneous == true`, then the analysis will also consider instantaneous interactions between
the variables.

If `maxlag < τcutoff` is an integer, `maxlag` is taken as the maximum allowed embedding lag. If `maxlag` is a float, 
then the maximum embedding lag is taken as `minimum(maximum([length.(source); length.(target); length.(cond)])*maxlag, τcutoff)`.
`τmin` is the minimum amount of variables to include in the candidate set (that is, each variable is embedded 
with lags `startlag:-1:-τmin`).
"""
function pa_candidate_variables(source, target, cond;
        k::Int = 1,
        include_instantaneous = true,
        method_delay = "ac_min",
        τmin::Int = 2,
        maxlag::Union{Int, Float64} = 0.05, 
        τcutoff::Int = 100)
    
    @assert τmin > 0

    # Generate candidate set. First, estimate a suitable maximum embedding lag for each 
    # of the input variables. We need at least two variables to allow exlusion of the `k`th
    # lag, so make sure we don't get only one variable from the delay method by using `τmin`.
    delay_τs = search_τs(source, target, cond, maxlag = maxlag, τcutoff = τcutoff)
    τsmax_source = [max(estimate_delay(s, method_delay, delay_τs), τmin) for s in source]
    τsmax_target = [max(estimate_delay(t, method_delay, delay_τs), τmin) for t in target]
    τsmax_cond = [max(estimate_delay(c, method_delay, delay_τs), τmin) for c in cond]

    # Next, compute the embeddings lags for each input variable, making sure that 
    # the prediction lag `k` does not overlap with the embedding lag.
    τs_source = [bbnue_embedding_lags(include_instantaneous, τ, k) for τ in τsmax_source]
    τs_target = [bbnue_embedding_lags(include_instantaneous, τ, k) for τ in τsmax_target]
    τs_cond = [bbnue_embedding_lags(include_instantaneous, τ, k) for τ in τsmax_cond]
    
    # Forward prediction lags.
    ks_T⁺ = [k for i in 1:length(target)]
    ks_T⁻ = [-k for i in 1:length(target)]
    js_T⁺ = [i for i in length(τs_source)+1:length(τs_source)+length(τs_target)]
    js_T⁻ = [i for i in length(τs_source)+1:length(τs_source)+length(τs_target)]

    τs_Ω = [τs_source..., τs_target..., τs_cond...]
    js_Ω = [[i for x in 1:length(τs[i])] for i = 1:length(τs)]

    return [τs_Ω...,], [js_Ω...,], ks_T⁺, ks_T⁻, js_T⁺, js_T⁻
end

function pa_candidate_variables(source, target; 
        k::Int = 1, 
        include_instantaneous = true,
        method_delay = "ac_min",
        τmin::Int = 2,
        maxlag::Union{Int, Float64} = 0.05, 
        τcutoff::Int = 100)
    
    # Generate candidate set. First, estimate a suitable maximum embedding lag for each 
    # of the input variables. We need at least two variables to allow exlusion of the `k`th
    # lag, so make sure we don't get only one variable from the delay method by using `τmin`.
    delay_τs = search_τs(source, target, cond, maxlag = maxlag, τcutoff = τcutoff)
    τsmax_source = [max(estimate_delay(s, method_delay, delay_τs), τmin) for s in source]
    τsmax_target = [max(estimate_delay(t, method_delay, delay_τs), τmin) for t in target]

    # Next, compute the embeddings lags for each input variable, making sure that 
    # the prediction lag `k` does not overlap with the embedding lag.
    τs_source = [bbnue_embedding_lags(include_instantaneous, τ, k) for τ in τsmax_source]
    τs_target = [bbnue_embedding_lags(include_instantaneous, τ, k) for τ in τsmax_target]
    
    ks_T⁺ = [k for i in 1:length(target)] 
    ks_T⁻ = [-k for i in 1:length(target)]
    js_T⁺ = [i for i in length(τs_source)+1:length(τs_source)+length(τs_target)]
    js_T⁻ = [i for i in length(τs_source)+1:length(τs_source)+length(τs_target)]

    τs_Ω = [τs_source..., τs_target...,]
    js_Ω = [[i for x in 1:length(τs[i])] for i = 1:length(τs)]

    return [τs_Ω...,], [js_Ω...,], ks_T⁺, ks_T⁻, js_T⁺, js_T⁻
end

function pa_embed_variables(source, target; 
        η::Int = 1, 
        include_instantaneous = true,
        method_delay = "mi_min",
        maxlag::Union{Int, Float64} = 0.05, 
        τcutoff::Int = 100)
    
    τs, js, ks_T⁺, ks_T⁻, js_T⁺, js_T⁻ = pa_candidate_variables(source, target, k = η, 
        τcutoff = τcutoff, maxlag = maxlag, method_delay = method_delay, 
        include_instantaneous = include_instantaneous)
    
    # TODO: This is more efficient if not using datasets. Re-do manually.
    data = Dataset([source..., target...,]...,)
    ℰ = genembed(data, ((τs...)...,), ((js...)...,))
    
    n_candidates = sum(length.(τs[1:end-2])) # the last two are target future and past
    n_T⁺ = length(τs[end-1]) 
    n_T⁻ = length(τs[end])

    nvars_total = length(((τs...)...,))

    # Separate time series for candidate variables (Ω), target future (T⁺), and target past (T⁻)
    idxs_candidates = 1:n_candidates
    idxs_target_future = n_candidates+1:nvars_total-n_T⁻
    idxs_target_past = nvars_total-n_T⁻+1:nvars_total
    Ω = [ℰ[:, i] for i = idxs_candidates]
    T⁺ = ℰ[:, idxs_target_future]
    T⁻ = ℰ[:, idxs_target_past]

    # @show τs
    # @show ((τs...)...,)
    # @show ((js...)...,)
    # @show n_candidates
    # @show idxs_candidates
    # @show idxs_target_future
    # @show idxs_target_past

    idxs_source = 1:length(source)
    idxs_target = length(source)+1:length(source)+length(target)
    idxs_cond = Int[]

    return Ω, T⁺, T⁻, τs, js, idxs_source, idxs_target, idxs_cond
end

"""
    minimize_conditional_entropy(T⁺, 𝒮, Ω) → idx

Minimize conditional entropy (CE) between the future of the target variable `T⁺`,
and (`𝒮`, Ωᵢ), where `Ωᵢ ∈ Ω` is the i-th remaining candidate variable.

Returns `i` such that `Ω[i]` is the variable that minimizes CE when combined with 
`𝒮`.
"""
function minimize_conditional_entropy(T⁺, 𝒮, Ω)
    n_remaining_candidates = length(Ω)
    ces = zeros(n_remaining_candidates)

    for i = 1:n_remaining_candidates
        if k == 1 || length(𝒮) == 0
            Cᵢ = Ω[i]
            CMI_T⁺_Cᵢ = 
                genentropy(Dataset(T⁺, Dataset(Cᵢ)), est, q = q, base = base) - 
                genentropy(Dataset(Cᵢ), est, q = q, base = base)
        else
            Cᵢ = [Ω[i], 𝒮...]
            CMI_T⁺_Cᵢ = 
                genentropy(Dataset(T⁺, Dataset(Cᵢ...,)), est, q = q, base = base) - 
                genentropy(Dataset(Cᵢ...,), est, q = q, base = base)
        end
        CMIs_between_T⁺_and_candidates[i] = CMI_T⁺_Cᵢ
    end

    return idx = findfirst(x -> x == minimum(ces), ces)
end



function optim_pa(Ω, T⁺, τs, js, idxs_source, idxs_target, idxs_cond, est; 
        uq = 0.95, nsurr = 100, q = 1, base = 2)
    
    τs_comb = [(τs...)...,]
    js_comb = [(js...)...,]
    
    npts = length(T⁺)
    n_candidate_variables = length(Ω)
    
    𝒮 = Vector{Vector{Float64}}(undef, 0)
    𝒮_τs = Vector{Int}(undef, 0)
    𝒮_js = Vector{Int}(undef, 0)
    
    k = 1
    while k <= n_candidate_variables
        # Find optimal variable for this `k`        
        idx = minimize_conditional_entropy(T⁺, 𝒮, Ω)
        Cₖ = Ω[idx]; sₖ = surrogenerator(Cₖ, CircShift(collect(1:npts - 1)))

        # Test significance of this candidate by using a random permutation test
        pa_surr = zeros(nsurr)
    
        if k == 1
            cmiₖ = CMIs_between_T⁺_and_candidates[idx]

            for i = 1:nsurr
                surr_Cₖ = sₖ() # Surrogate version of Cₖ
                CMI_permutations[i] = mutualinfo(T⁺, surr_Cₖ, est)
            end
        else
            # Precompute terms that do not change during surrogate loop
            H_T⁺_𝒮 = genentropy(Dataset(T⁺, Dataset(𝒮...,)), est, q = q, base = base)
            
            # ORIGIANL TE
            H_𝒮 = genentropy(Dataset(𝒮...), est, q = q, base = base)
            cmiₖ = H_T⁺_𝒮 + 
                    genentropy(Dataset([Cₖ, 𝒮...,]...,), est, q = q, base = base) - 
                    genentropy(Dataset(T⁺, Dataset([Cₖ, 𝒮...,]...,)), est, q = q, base = base) - 
                    H_𝒮

            for i = 1:nsurr
                surr_Cₖ = s() # Surrogate version of Cₖ
                CMI_permutations[i] = H_T⁺_𝒮 + 
                    genentropy(Dataset([surr_Cₖ, 𝒮...]...,), est, q = q, base = base) - 
                    genentropy(Dataset(T⁺, Dataset([surr_Cₖ, 𝒮...]...,)), est, q = q, base = base) - 
                    H_𝒮
            end
            
        end
       # If the candidate passes the significance test
        if cmiₖ > quantile(CMI_permutations, uq)
            # Add the candidate to list of selected candidates
            push!(𝒮, Cₖ)
            push!(𝒮_τs, τs_comb[idx])
            push!(𝒮_js, js_comb[idx])
            
            # Delete the candidate from the list of remaining candidates
            deleteat!(Ω, idx)
            deleteat!(τs_comb, idx)
            deleteat!(js_comb, idx)

            k = k + 1
        else 
            k = n_candidate_variables + 1
        end
    end
    
    
    # No variables were selected
    if length(𝒮) == 0
        return 0.0, Int[], Int[], idxs_source, idxs_target, idxs_cond
    end
    
    # No variables were selected from the source process
    n_source_vars_picked = count(x -> x ∈ idxs_source, 𝒮_js)
    if n_source_vars_picked == 0
        return 0.0, Int[], Int[], idxs_source, idxs_target, idxs_cond
    end
    
    # No variables were selected from the target or conditional processes.
    𝒮_nonX = [ts for (ts, j) in zip(𝒮, 𝒮_js) if j ∉ idxs_source]
    if length(𝒮_nonX) == 0
        return 0.0, Int[], Int[], idxs_source, idxs_target, idxs_cond
    end
        
    CE2 = genentropy(Dataset(T⁺, Dataset(𝒮...,)), est, base = base, q = q) - 
        genentropy(Dataset(𝒮...,), est, base = base, q = q)
    
    CE1 = genentropy(Dataset(T⁺, Dataset(𝒮_nonX...,)), est, base = base, q = q) - 
        genentropy(Dataset(𝒮_nonX...,), est, base = base, q = q)
    
    CMI = CE1 - CE2
    return CMI, 𝒮_js, 𝒮_τs, idxs_source, idxs_target, idxs_cond
    
end

# TODO: verify selection in separate method?