using DelayEmbeddings, Entropies

export preprocess_pa
function preprocess_pa(targets, sources; τmax = 2, ηmax = 5)
    Lt = length(targets)
    Ls = length(sources)
    
    τs = 0:-1:-τmax
    ηs = 1:ηmax
    #@show τs |> collect
    #@show ηs |> collect


    # js for `genembed` are ordered as follows:
    # target1 = var1 
    # target2 = var2
    # ...
    # target_n = var_n
    # source1 = var_n+1
    # source2 = var_n+2
    # ....

    js = [[repeat([i], length(ηs)) for (i, t) in enumerate(targets)]...; 
            [repeat([i], length(τs)) for (i, t) in enumerate(targets)]...; 
            [repeat([i+length(targets)], length(τs)) for (i, t) in enumerate(sources)]...]

    # The embedding lags corresponding to `js`
    lags = [[ηs for t in targets]...; 
            [τs for t in targets]...; 
            [τs for s in sources]...]

    allvars = Dataset(targets..., sources...,)
    D = genembed(allvars, lags, js)

    # Select variables
    nη = length(ηs)
    nτ = length(τs)
    
    # Select prediction variables
    ts_pred = [D[:, i] for i = 1:size(D, 2)]
    ts_targets = [[ts_pred[i + j*nη] for j = 0:Lt-1] for i = 1:nη]
    
    # Select candidate variables
    ts_cond = [D[:, i] for i = nη*Lt+1:size(D, 2)]
    idxs_candidates = [i:i+(nτ-1) for i in [k + (nτ - 1)*(k - 1) for k in 1:(Ls+Lt)]]
    ts_candidates = [ts_cond[idxs] for idxs in idxs_candidates]

    @show js   |> collect
    @show lags |> collect
    ts_targets, ts_candidates
end

export te_aut

# Assumes pre-embedded time series

"""
    te_aut(targets, candidates, est, q = 1.0, base = 2)

An implementation of the bootstrap-based non-uniform embedding (BB-NUE) algorithm for transfer entropy estimation.

Takes a set of `targets` and `candidates` and estimates transfer entropy from `candidates` to `targets`, selecting 
only the variables from `candidates` that are the most relevant and least redundant past variables.

Assumes that variables have been pre-processed using `preprocess_pa`, which returns two arrays: 
- `targets`, where `targets[η]` contains a vector of variables corresponding to the target at predition lag `η` (may be multidimensional).
- `candidates`, where `candidates[k]` contains a vector of lagged variables for the `k`-th candidate variables. There are `length(τs)` variables
    for each input variable.
    
[Shahsavari2020]: Shahsavari Baboukani, P., Graversen, C., Alickovic, E., & Østergaard, J. (2020). Estimating Conditional Transfer Entropy in Time Series Using Mutual Information and Nonlinear Prediction. Entropy, 22(10), 1124. https://www.mdpi.com/1099-4300/22/10/1124/htm
"""
function te_aut(targets, candidates, est, q = 1.0, base = 2)
    # Infer prediction lags from the number of target variables (targets is always a vector of vectors)
    ηs = 1:length(targets)

    # Infer lags from first set of candidates (candidates is always a vector of vectors)
    τs = 0:-1:-length(candidates)[1] - 1    
    #@show τs |> collect, length(τs)
    #@show ηs |> collect, length(ηs)

    candidates_flattened = Iterators.flatten(candidates)

    N = length(candidates)

    # For every η, find optimal variables
    for η in ηs
        # Already chosen variables.
        θs = []

        k = 1; kmax = 5
        while k < kmax
            #target = Dataset(θs...)
            # The first chosen variable is the target variable.
            #push!(θs, targets[η]...)

            Ntot = sum(length.(candidates))
            Hs = zeros(Ntot)
            Hs_idxs_flattened = zeros(Int, Ntot)
            Hs_idxs = zeros(Int, Ntot)
            Hs_lags = zeros(Int, Ntot)

            # For the i = 1, 2..., 
            ct = 1
            for i = 1:N
                # The candidates for the i-th candidate
                Cᵢ = candidates[i]
                Hsᵢ = zeros(length(Cᵢ))
                
                for k = 1:length(Cᵢ)
                    # H(θ | Cᵢ[k]) = H(θ, Cᵢ[k]) - H(Cᵢ[k]) (chain rule for entropy)
                    Hθ_and_Cᵢk = genentropy(Dataset(targets[η]..., Cᵢ[k]), est, q = q, base = base)
                    HCᵢk = genentropy(Dataset(Cᵢ[k]), est, q = q, base = base)
                    Hθ_given_Cᵢk = Hθ_and_Cᵢk - HCᵢk
                    Hsᵢ[k] = Hθ_given_Cᵢk 
                    Hs[ct] = Hθ_given_Cᵢk
                    Hs_idxs[ct] = i
                    Hs_idxs_flattened[ct] = ct
                    Hs_lags[ct] = τs[k]

                    ct += 1
                end
            end

            # Identify the candidate variable which reduces conditional entropy the most (maximizes the 
            # shared information between `θs` and the candidate)
            #@show Hs
            #@show Hs_idxs
            #@show Hs_lags

            idx = findmin(Hs)[2]
            j = Hs_idxs[idx]
            varidx = Hs_idxs_flattened[idx]
            τidx = -Hs_lags[idx] + 1
            chosen_candidate = candidates[j][τidx]


            println("The most informative candidate was candidate #$j at lag $(Hs_lags[idx])")


            passes = rand([true, false])

            # if passes
            #     deleteat!(candidates[j], τidx)
            #     println("\t The candidate passed the significance test. Removes its status as candidate. $Ntot candidates remaining.")
            # else
            #     println("\t The candidate failed the significance test")
            # end
            cmi = 0.0

            if k == 1 # No variables have been selected (θ == ∅), so CMI reduces to MI
                candidate = Dataset(chosen_candidate)
                target = Dataset(targets[η]...)
                Hc = genentropy(Dataset(candidate), est; base = base, q = q)
                Ht = genentropy(Dataset(target), est; base = base, q = q)
                Hct = genentropy(Dataset(candidate, target), est; base = base, q = q)
                cmi = Hc + Ht - Hct
            else # At least one variable has been selected, so we can condition on θ

            end
            @show cmi
            # Make sure we terminate at some point
            k += 1
        end

        println("Done")
        

        
    end

end