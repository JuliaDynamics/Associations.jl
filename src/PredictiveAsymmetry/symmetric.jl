

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