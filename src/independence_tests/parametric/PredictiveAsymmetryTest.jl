using HypothesisTests

"""
    PredictiveAsymmetryTest <: IndependenceTest
    PredictiveAsymmetryTest(measure::PredictiveAsymmetryTest, est; f = 1.0)

An independence test based on the [`PredictiveAsymmetry`](@ref) directional association
measure. When used with [`independence`](@ref), a p-value is returned.
"""
struct PredictiveAsymmetryTest{M, E, F} <: IndependenceTest{M}
    measure::M
    est::E
    f::F

    function PredictiveAsymmetryTest(measure::M = PredictiveAsymmetry(), est::E = Zhu1();
            f::F = 1) where {M, E, F}
        return new{M, E, F}(measure, est, f)
    end
end

struct PredictiveAsymmetryTestResult{A, TF, TB, T, P} <: IndependenceTestResult
    n_vars::Int # 2 vars = pairwise, 3 vars = conditional
    ΔA::A
    tes_fw::TF
    tes_bw::TB
    μ0::T
    pvalue::P
end

pvalue(r::PredictiveAsymmetryTestResult) = r.pvalue

function Base.show(io::IO, test::PredictiveAsymmetryTestResult)
    print(io,
        """\
        `PredictiveAsymmetryTestResult` independence test
        $(null_hypothesis_text(test))
        $(pvalue_text_summary(test))
        """
        )
end

function independence(test::PredictiveAsymmetryTest, source, target)
    (; ηs, normalize, f, base, dTf, dT, dS, dC, τT, τS, τC) = test.measure
    est = test.est
    check_ηs(ηs)
    Nη = length(ηs)

    te_fws = zeros(Nη)
    te_bws = zeros(Nη)
    ΔA = zeros(Nη)

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
            ΔA[i] = 𝔸ᵢ / (f*avg_te)
        else
            𝔸ᵢ = (sum(te_fws[1:i]) - sum(te_bws[1:i])) / η
            ΔA[i] = 𝔸ᵢ
        end
    end

    # Use mean transfer entropy as threshold for significance
    μ0 = test.f * mean([te_fws; te_bws])
    p = HypothesisTests.pvalue(OneSampleTTest(ΔA, μ0), tail = :right)

    return PredictiveAsymmetryTestResult(2, ΔA, te_fws, te_bws, μ0, p)
end

function independence(test::PredictiveAsymmetryTest, source, target, cond)
    (; ηs, normalize, f, base, dTf, dT, dS, dC, τT, τS, τC) = test.measure
    est = test.est
    check_ηs(ηs)
    Nη = length(ηs)

    te_fws = zeros(Nη)
    te_bws = zeros(Nη)
    ΔA = zeros(Nη)

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
            ΔA[i] = 𝔸ᵢ / (f*avg_te)
        else
            𝔸ᵢ = (sum(te_fws[1:i]) - sum(te_bws[1:i])) / η
            ΔA[i] = 𝔸ᵢ
        end
    end

   # Use mean transfer entropy as threshold for significance
   μ0 = test.f * mean([te_fws; te_bws])
   p = HypothesisTests.pvalue(OneSampleTTest(ΔA, μ0), tail = :right)

   return PredictiveAsymmetryTestResult(3, ΔA, te_fws, te_bws, μ0, p)
end
