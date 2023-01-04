import ComplexityMeasures: symbolize_for_dispersion

function marginal_encodings(est::SymbolicPermutation{m}, x::AbstractVector...) where {m}
    return [encode.(Ref(est.encoding), embed(xᵢ, m, est.τ).data) for xᵢ in x]
end

function marginal_encodings(est::Dispersion{m}, x::AbstractVector...) where {m}
    return [symbolize_for_dispersion(est, xᵢ) for xᵢ in x]
end
