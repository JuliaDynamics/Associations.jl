

export pa_naive
# PA = H(Ȳ⁻, Y⁺, X⁻) - H(Ȳ⁻, Y⁺) - H(Ȳ⁺, Y⁻, X⁻) + H(Ȳ⁺, Y⁻)
function pa_naive(source, target, est::VisitationFrequency{RectangularBinning{Int}}, ηs;
        d𝒯 = 1, dT = 1, dS = 1, τT = -1, τS = -1, base = 2, q = 1)
    @assert τT < 0
    @assert τS < 0 
    τs⁻ = 0:-1:τS*(dS - 1)
    νs⁻ = 0:-1:τT*(dT - 1)
    νs⁺ = .-(0:-1:τT*(dT - 1))
    data = Dataset(source, target)
    lags = [τs⁻..., νs⁻..., νs⁺..., ηs..., .-(ηs)...,]
    js = [repeat([1], length(τs⁻))...,repeat([2], length(νs⁻)*2)..., repeat([2], length(ηs)*2)..., ]
    E = genembed(data, lags, js)

    nη = length(ηs)
    X⁻ = Dataset(E[:, 1:dS])
    Y⁻ = Dataset(E[:, dS+1:dS+dT])
    Y⁺ = Dataset(E[:, dS+dT+1:dS+dT+dT])
    n = dimension(X⁻)+dimension(Y⁻)+dimension(Y⁺)
    Ȳ⁻ = Dataset(E[:, n+1:n+nη])
    Ȳ⁺ = Dataset(E[:, n+nη+1:end])
    
    pas = zeros(nη)
    #@show 1:dS, dS+1:dS+dT, dS+dT+1:dS+dT+dT, n+1:n+length(ηs), n+length(ηs)+1:length(lags)
    # Ȳ⁺Y⁻X⁻s = zeros(nη)
    # Ȳ⁻Y⁺X⁻s = zeros(nη)
    # Ȳ⁺Y⁻s = zeros(nη)
    # Ȳ⁻Y⁺s = zeros(nη)

    for η in ηs
        Ȳ⁻Y⁺ = Dataset(Dataset(Ȳ⁻[:, η]), Y⁺)
        Ȳ⁻Y⁺X⁻ = Dataset(Ȳ⁻Y⁺, X⁻)

        Ȳ⁺Y⁻ = Dataset(Dataset(Ȳ⁺[:, η]), Y⁻)
        Ȳ⁺Y⁻X⁻ = Dataset(Ȳ⁺Y⁻, X⁻)

        nbins = est.binning.ϵ

        # Ȳ⁺Y⁻X⁻s[η] = genentropy(Ȳ⁺Y⁻X⁻, est, base = base, q = q)
        # Ȳ⁻Y⁺X⁻s[η] = genentropy(Ȳ⁻Y⁺X⁻, est, base = base, q = q)
        # Ȳ⁺Y⁻s[η] = genentropy(Ȳ⁺Y⁻, est, base = base, q = q) 
        # Ȳ⁻Y⁺s[η] = genentropy(Ȳ⁻Y⁺, est, base = base, q = q)
        
        pas[η] = 
            genentropy(Ȳ⁻Y⁺X⁻, est, base = base, q = q) -
            genentropy(Ȳ⁻Y⁺, est, base = base, q = q)  -
            genentropy(Ȳ⁺Y⁻X⁻, est, base = base, q = q) +
            genentropy(Ȳ⁺Y⁻, est, base = base, q = q)  
            
    end

    # Maximal entropy of a system with N states is log(N)
    return pas#, Ȳ⁺Y⁻X⁻s, Ȳ⁻Y⁺X⁻s, Ȳ⁺Y⁻s, Ȳ⁻Y⁺s
end

export pa_crapload_inverse
function pa_crapload_inverse(source, target, est::VisitationFrequency{RectangularBinning{Int}}, ηs;
        d𝒯 = 1, dT = 1, dS = 1, τT = -1, τS = -1, base = 2, q = 1)
    @assert τT < 0
    @assert τS < 0 
    τs⁻ = 0:-1:τS*(dS - 1)
    νs⁻ = 0:-1:τT*(dT - 1)
    νs⁺ = .-(0:-1:τT*(dT - 1))
    data = Dataset(source, target)
    lags = [τs⁻..., νs⁻..., νs⁺..., ηs..., .-(ηs)...,]
    js = [repeat([1], length(τs⁻))...,repeat([2], length(νs⁻)*2)..., repeat([2], length(ηs)*2)..., ]
    E = genembed(data, lags, js)

    nη = length(ηs)
    X⁻ = Dataset(E[:, 1:dS])
    Y⁻ = Dataset(E[:, dS+1:dS+dT])
    Y⁺ = Dataset(E[:, dS+dT+1:dS+dT+dT])
    n = dimension(X⁻)+dimension(Y⁻)+dimension(Y⁺)
    Ȳ⁻ = Dataset(E[:, n+1:n+nη])
    Ȳ⁺ = Dataset(E[:, n+nη+1:end])
    
    pas = zeros(nη)
    #@show 1:dS, dS+1:dS+dT, dS+dT+1:dS+dT+dT, n+1:n+length(ηs), n+length(ηs)+1:length(lags)

    for η in ηs
        Y⁺X⁻ = Dataset(Y⁺, X⁻)
        Y⁻X⁻ = Dataset(Y⁻, X⁻)

        nbins = est.binning.ϵ
        pas[η] = 
            genentropy(Y⁺X⁻, est, base = base, q = q) -
            genentropy(Y⁻X⁻, est, base = base, q = q)
    end

    # Maximal entropy of a system with N states is log(N)
    return pas
end



# PA = H(Ȳ⁻, Y⁺, X⁻) - H(Ȳ⁻, Y⁺) - H(Ȳ⁺, Y⁻, X⁻) + H(Ȳ⁺, Y⁻)
export pa_te_based
function pa_te_based(source, target, est::VisitationFrequency{RectangularBinning{Int}}, ηs;
        d𝒯 = 1, dT = 1, dS = 1, τT = -1, τS = -1, base = 2, q = 1)
    @assert τT < 0
    @assert τS < 0 
    τs⁻ = 0:-1:τS*(dS - 1)
    νs⁻ = 0:-1:τT*(dT - 1)
    νs⁺ = .-(0:-1:τT*(dT - 1))
    data = Dataset(source, target)
    lags = [τs⁻..., νs⁻..., νs⁺..., ηs..., .-(ηs)...,]
    js = [repeat([1], length(τs⁻))...,repeat([2], length(νs⁻)*2)..., repeat([2], length(ηs)*2)..., ]
    E = genembed(data, lags, js)

    nη = length(ηs)
    X⁻ = Dataset(E[:, 1:dS])
    Y⁻ = Dataset(E[:, dS+1:dS+dT])
    Y⁺ = Dataset(E[:, dS+dT+1:dS+dT+dT])
    n = dimension(X⁻)+dimension(Y⁻)+dimension(Y⁺)
    Ȳ⁻ = Dataset(E[:, n+1:n+nη])
    Ȳ⁺ = Dataset(E[:, n+nη+1:end])
    
    pas = zeros(nη)

    for η in ηs
        Ȳ⁺Y⁻ = Dataset(Dataset(Ȳ⁺[:, η]), Y⁻)
        X⁻Y⁻ = Dataset(X⁻, Y⁻)
        Ȳ⁺Y⁻X⁻ = Dataset(Ȳ⁺Y⁻, X⁻)

        Ȳ⁻Y⁺ = Dataset(Dataset(Ȳ⁻[:, η]), Y⁺)
        X⁻Y⁺ = Dataset(X⁻, Y⁺)
        Ȳ⁻Y⁺X⁻ = Dataset(Ȳ⁻Y⁺, X⁻)
        
        nbins = est.binning.ϵ
        # When starting from TE and flipping embedding
        # H(Ȳ⁺, Y⁻) + H(X⁻, Y⁻) - H(Y⁻) - H(X⁻, Y⁻, Ȳ⁺) - 
        #   H(Ȳ⁻, Y⁺) - H(X⁻, Y⁺) + H(Y⁺) + H(X⁻, Y⁺, Ȳ⁻) - 
        pas[η] = 
            genentropy(Ȳ⁺Y⁻, est, base = base, q = q) + #/ (nbins * dimension(Ȳ⁺Y⁻)) + 
            genentropy(X⁻Y⁻, est, base = base, q = q) - #/ (nbins * dimension(X⁻Y⁻)) -
            genentropy(Y⁻, est, base = base, q = q) - #/ (nbins * dimension(Y⁻)) -
            genentropy(Ȳ⁺Y⁻X⁻, est, base = base, q = q) - #/ (nbins * dimension(Ȳ⁺Y⁻X⁻)) - 
            genentropy(Ȳ⁻Y⁺, est, base = base, q = q) - #/ (nbins * dimension(Ȳ⁻Y⁺)) - 
            genentropy(X⁻Y⁺, est, base = base, q = q) + #/ (nbins * dimension(X⁻Y⁺)) + 
            genentropy(Y⁺, est, base = base, q = q) + #/ (nbins * dimension(Y⁺)) +
            genentropy(Ȳ⁻Y⁺X⁻, est, base = base, q = q) #/ (nbins * dimension(Ȳ⁻Y⁺X⁻)) 
    end

    # Maximal entropy of a system with N states is log(N)
    return pas
end

struct JointEntropies
    Ȳ⁻Y⁺s
    Ȳ⁺Y⁻s
    X⁻Ȳ⁺Y⁻s
    X⁻Ȳ⁻Y⁺s
    X⁺Ȳ⁺Y⁻s
    X⁺Ȳ⁻Y⁺s
end

export relevant_joint_entropies
function relevant_joint_entropies(source, target, est::VisitationFrequency{RectangularBinning{Int}}, ηs;
        dT = 1, dS = 1, τT = -1, τS = -1, base = 2, q = 1)
    @assert τT < 0
    @assert τS < 0 

    τs⁻ = 0:τS:τS*(dS - 1)
    τs⁺ = .-(τs⁻)
    νs⁻ = 0:τT:τT*(dT - 1)
    νs⁺ = .-(νs⁻)
    nη = length(ηs)

    data = Dataset(source, target)
    lags = [τs⁻..., τs⁺..., νs⁻..., νs⁺..., ηs..., .-(ηs)...,]
    js = [repeat([1], dS*2)..., repeat([2], dT*2)..., repeat([2], nη*2)..., ]
    E = genembed(data, lags, js)

    X⁻ = Dataset(E[:, 1:dS])
    X⁺ = Dataset(E[:, dS+1:dS*2])
    Y⁻ = Dataset(E[:, dS*2+1:dS*2+dT])
    Y⁺ = Dataset(E[:, dS*2+dT+1:dS*2+dT*2])
    Ȳ⁻ = Dataset(E[:, dS*2+dT*2+1:dS*2+dT*2+nη])
    Ȳ⁺ = Dataset(E[:, dS*2+dT*2+nη+1:end])

    Ȳ⁻Y⁺s = zeros(nη)
    Ȳ⁺Y⁻s = zeros(nη)
    X⁻Ȳ⁺Y⁻s = zeros(nη)
    X⁻Ȳ⁻Y⁺s = zeros(nη)
    X⁺Ȳ⁺Y⁻s = zeros(nη)
    X⁺Ȳ⁻Y⁺s = zeros(nη)

    for η in ηs
        Ȳ⁻Y⁺ = Dataset(Dataset(Ȳ⁻[:, η]), Y⁺)
        Ȳ⁺Y⁻ = Dataset(Dataset(Ȳ⁺[:, η]), Y⁻)
        X⁻Ȳ⁺Y⁻ = Dataset(X⁻, Ȳ⁺Y⁻)
        X⁻Ȳ⁻Y⁺ = Dataset(X⁻, Ȳ⁻Y⁺)
        X⁺Ȳ⁺Y⁻ = Dataset(X⁺, Ȳ⁺Y⁻)
        X⁺Ȳ⁻Y⁺ = Dataset(X⁺, Ȳ⁻Y⁺)

        Ȳ⁻Y⁺s[η] = genentropy(Ȳ⁻Y⁺, est, base = base, q = q)
        Ȳ⁺Y⁻s[η] = genentropy(Ȳ⁺Y⁻, est, base = base, q = q)
        X⁻Ȳ⁺Y⁻s[η] = genentropy(X⁻Ȳ⁺Y⁻, est, base = base, q = q)
        X⁻Ȳ⁻Y⁺s[η] = genentropy(X⁻Ȳ⁻Y⁺, est, base = base, q = q)
        X⁺Ȳ⁺Y⁻s[η] = genentropy(X⁺Ȳ⁺Y⁻, est, base = base, q = q)
        X⁺Ȳ⁻Y⁺s[η] = genentropy(X⁺Ȳ⁻Y⁺, est, base = base, q = q)
    end
    #@show lags
    #@show js
    #@show "X⁻", 1:dS
    #@show "X⁺", dS+1:dS*2
    #@show "Y⁻", dS*2+1:dS*2+dT
    #@show "Y⁺", dS*2+dT+1:dS*2+dT*2
    #@show "Ȳ⁻", dS*2+dT*2+1:dS*2+dT*2+nη
    #@show "Ȳ⁺", dS*2+dT*2+nη+1:length(lags)
    #eturn [js lags
    return JointEntropies(Ȳ⁻Y⁺s, Ȳ⁺Y⁻s, X⁻Ȳ⁺Y⁻s, X⁻Ȳ⁻Y⁺s, X⁺Ȳ⁺Y⁻s, X⁺Ȳ⁻Y⁺s)
end


struct CondEntropies
    X⁻_Ȳ⁺Y⁻
    X⁻_Ȳ⁻Y⁺
    X⁺_Ȳ⁺Y⁻
    X⁺_Ȳ⁻Y⁺
end

export relevant_cond_entropies
function relevant_cond_entropies(source, target, est::VisitationFrequency{RectangularBinning{Int}}, ηs;
        dT = 1, dS = 1, τT = -1, τS = -1, base = 2, q = 1)
    @assert τT < 0
    @assert τS < 0 
    H = relevant_joint_entropies(source, target, est, ηs, 
        dT = dT, dS = dS, τT = τT, τS = τS, base = base, q = q)
    
    CondEntropies(
        H.X⁻Ȳ⁻Y⁺s .- H.Ȳ⁻Y⁺s,
        H.X⁻Ȳ⁻Y⁺s .- H.Ȳ⁻Y⁺s,
        H.X⁺Ȳ⁺Y⁻s .- H.Ȳ⁺Y⁻s,
        H.X⁺Ȳ⁻Y⁺s .- H.Ȳ⁻Y⁺s
    )
end

function A₁(source, target, est::VisitationFrequency{RectangularBinning{Int}}, ηs;
        dT = 1, dS = 1, τT = -1, τS = -1, base = 2, q = 1)
    H = relevant_cond_entropies(source, target, est, ηs, 
        dT = dT, dS = dS, τT = τT, τS = τS, base = base, q = q)

    H.X⁺_Ȳ⁺Y⁻ .- H.X⁻_Ȳ⁺Y⁻
end

function A₂(source, target, est::VisitationFrequency{RectangularBinning{Int}}, ηs;
    dT = 1, dS = 1, τT = -1, τS = -1, base = 2, q = 1)
    H = relevant_cond_entropies(source, target, est, ηs, 
        dT = dT, dS = dS, τT = τT, τS = τS, base = base, q = q)

    H.X⁺_Ȳ⁻Y⁺ .- H.X⁻_Ȳ⁻Y⁺
end