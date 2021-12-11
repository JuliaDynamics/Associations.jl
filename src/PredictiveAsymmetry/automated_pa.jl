using Statistics

function pa_bbnue(x, y, est, ηs; 
        nsurr = 100, uq = 0.95, 
        include_instantaneous = true, 
        method_delay = "ac_min", 
        maxlag::Union{Int, Float64} = 0.05)
    
    tes_fw = zeros(length(ηs))
    tes_bw = zeros(length(ηs))

    for (i, η) in enumerate(ηs)
        fw_te, fw_js, fw_τs = bbnue(x, y, est, η = η;
            nsurr = nsurr, uq = uq, 
            include_instantaneous = include_instantaneous, 
            method_delay = method_delay, maxlag = maxlag, q = q)

        bw_te, bw_js, bw_τs = bbnue(x, y, est, η = -η; 
            nsurr = nsurr, uq = uq, 
            include_instantaneous = include_instantaneous, 
            method_delay = method_delay, maxlag = maxlag, q = q)
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

    mean(opt_pa_bbnue(x, y, est, ηs; 
        nsurr = nsurr, uq = uq, 
        include_instantaneous = include_instantaneous,
        method_delay = method_delay, maxlag = maxlag))
end