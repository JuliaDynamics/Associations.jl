

@recipe function f{CT <: AbstractPredictiveAsymmetryTest, DT, RT}(analysis::CausalAnalysis{CT,DT,RT})
    pas = analysis.result
    ymax = maximum(pas) * 1.1
    ηs = get_ηs(analysis.test)
    ηs = ηs[ηs .> 0]
    
    if CT isa NormalisedPredictiveAsymmetryTest
        ylabel --> L"$\mathbb{A}(\eta)$"
    elseif CT isa PredictiveAsymmetryTest
        ylabel --> L"$\mathcal{A}(\eta)$"
    end
    seriestype --> :line
    xlabel --> L"$\eta$"
    ylabel --> L"$\mathbb{A}(\eta)$"
    label --> ""
    size --> (382, 300)
    fg_legend --> :transparent
    bg_legend --> :transparent
       
    (ηs, pas)
end
