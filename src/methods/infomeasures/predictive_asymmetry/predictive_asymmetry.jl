export PredictiveAsymmetry
export PredictiveAsymmetryTest

"""
    PredictiveAsymmetry <: AssociationMeasure
    PredictiveAsymmetry(ηs = 1:15; normalize = false, f = 1.0,
        dTf = 1, dT = 1, dS = 1, τT = -1, τS = -1, base = 2)

The predictive asymmetry measure [Haaga2020](@cite).

!!! info "Experimental!"
    This is a method that does not yet appear in a peer-reviewed scientific journal.
    Feel free to use, but consider it experimental for now.
"""
Base.@kwdef struct PredictiveAsymmetry{E, B, F, D1, D2, D3, D4, T1, T2, T3} <: AssociationMeasure
    ηs::E = 1:15
    normalize::Bool = false
    f::F = 1.0
    base::B = 2
    dTf::D1 = 1
    dT::D2 = 1
    dS::D3 = 1
    dC::D4 = 1
    τT::T1 = -1
    τS::T2 = -1
    τC::T3 = -1
end

max_inputs_vars(::PredictiveAsymmetry) = 3

const PA_ESTIMATORS = Union{
    ProbabilitiesEstimator,
    DifferentialEntropyEstimator,
    MutualInformationEstimator,
    ConditionalMutualInformationEstimator,
    TransferEntropyEstimator
    }

function association(measure::PredictiveAsymmetry, est::PA_ESTIMATORS, source, target)
    (; ηs, normalize, f, base, dTf, dT, dS, dC, τT, τS, τC) = measure

    check_ηs(ηs)
    Nη = length(ηs)

    te_fws = zeros(Nη)
    te_bws = zeros(Nη)
    𝔸s = zeros(Nη)

    for (i, η) in enumerate(ηs)
        emb_fw = EmbeddingTE(dTf = dTf, dT = dT, dS = dS, τT = τT, τS = τS, ηTf = η)
        emb_bw = EmbeddingTE(dTf = dTf, dT = dT, dS = dS, τT = τT, τS = τS, ηTf = -η)
        te_fw = TEShannon(; base, embedding = emb_fw)
        te_bw = TEShannon(; base, embedding = emb_bw)
        te_fws[i] = transferentropy(te_fw, est, source, target)
        te_bws[i] = transferentropy(te_bw, est, source, target)

        if normalize
            𝔸ᵢ = (sum(te_fws[1:i]) - sum(te_bws[1:i])) / η
            avg_te = (sum(te_fws[1:i]) + sum(te_bws[1:i])) / (2*η)
            𝔸s[i] = 𝔸ᵢ / (f*avg_te)
        else
            𝔸ᵢ = (sum(te_fws[1:i]) - sum(te_bws[1:i])) / η
            𝔸s[i] = 𝔸ᵢ
        end
    end
    return 𝔸s
end

function association(measure::PredictiveAsymmetry, est::PA_ESTIMATORS, source, target, cond)
    (; ηs, normalize, f, base, dTf, dT, dS, dC, τT, τS, τC) = measure

    check_ηs(ηs)
    Nη = length(ηs)

    te_fws = zeros(Nη)
    te_bws = zeros(Nη)
    𝔸s = zeros(Nη)

    for (i, η) in enumerate(ηs)
        emb_fw = EmbeddingTE(dTf = dTf, dT = dT, dS = dS, dC = dC, τT = τT, τS = τS, τC = τC, ηTf = η)
        emb_bw = EmbeddingTE(dTf = dTf, dT = dT, dS = dS, dC = dC, τT = τT, τS = τS, τC = τC, ηTf = -η)
        te_fw = TEShannon(; base, embedding = emb_fw)
        te_bw = TEShannon(; base, embedding = emb_bw)
        te_fws[i] = transferentropy(te_fw, est, source, target, cond)
        te_bws[i] = transferentropy(te_bw, est, source, target, cond)

        if normalize
            𝔸ᵢ = (sum(te_fws[1:i]) - sum(te_bws[1:i])) / η
            avg_te = (sum(te_fws[1:i]) + sum(te_bws[1:i])) / (2*η)
            𝔸s[i] = 𝔸ᵢ / (f*avg_te)
        else
            𝔸ᵢ = (sum(te_fws[1:i]) - sum(te_bws[1:i])) / η
            𝔸s[i] = 𝔸ᵢ
        end
    end
    return 𝔸s
end
