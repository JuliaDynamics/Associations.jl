using Statistics
include("autoutils.jl")

export pa_bbnue, pa_mean_bbnue

function pa_bbnue(x, y, est, ηs; 
        nsurr = 100, uq = 0.95, 
        include_instantaneous = true, 
        method_delay = "ac_min", 
        maxlag::Union{Int, Float64} = 0.05,
        q = 1, base = 2)
    
    tes_fw = zeros(length(ηs))
    tes_bw = zeros(length(ηs))

    for (i, η) in enumerate(ηs)
        fw_te, fw_js, fw_τs = _pa_nue(x, y, est, η = η;
            nsurr = nsurr, uq = uq, 
            include_instantaneous = include_instantaneous, 
            method_delay = method_delay, maxlag = maxlag, 
            q = q, base = base)

        bw_te, bw_js, bw_τs = _pa_nue(x, y, est, η = -η; 
            nsurr = nsurr, uq = uq, 
            include_instantaneous = include_instantaneous, 
            method_delay = method_delay, maxlag = maxlag, 
            q = q, base = base)
        tes_fw[i] = fw_te
        tes_bw[i] = bw_te
    end
    
    pas_xy = [sum(tes_fw[1:i] .- tes_bw[1:i])/i for i = 1:length(ηs)]
end

function pa_mean_bbnue(x, y, est, ηs; 
        nsurr = 100, uq = 0.95, 
        include_instantaneous = true, 
        method_delay = "ac_min", 
        maxlag::Union{Int, Float64} = 0.05)

    mean(pa_bbnue(x, y, est, ηs; 
        nsurr = nsurr, uq = uq, 
        include_instantaneous = include_instantaneous,
        method_delay = method_delay, maxlag = maxlag))
end


function _pa_nue(source, target, est; 
    q = 0.95, base = 2,
    η = 1, nsurr = 100, uq = 0.95, 
    include_instantaneous = true, 
    method_delay = "ac_min", 
    maxlag::Union{Int, Float64} = 0.05)

Ω, Y⁺, τs, js, idxs_source, idxs_target, idxs_cond = 
    embed_candidate_variables(
        process_input(source), 
        process_input(target), 
        η = η, τexclude = η)

return optim_te(Ω, Y⁺, τs, js, idxs_source, idxs_target, idxs_cond, est, 
    q = q, base = base,
    nsurr = nsurr, uq = uq)
end