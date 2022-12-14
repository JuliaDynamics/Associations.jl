export RelativeEntropyShannonDifferential
export ShannonDivergenceDifferential

"""
    KLDivergenceDifferential  <: DivergenceDefinition

A directive to compute the differential Kullback-Leibler divergence (or Shannon relative
entropy).
"""
struct ShannonDivergenceDifferential  <: DivergenceDefinition end

"""
    RelativeEntropyShannonDifferential <: Divergence
    RelativeEntropyShannonDifferential(base = 2,
        definition::D = KLDivergenceDifferential())

[^Bulinski2021]:
    Bulinski, A., & Dimitrov, D. (2021). Statistical estimation of the Kullback-Leibler
    divergence. Mathematics, 9(5), 544.
"""
struct RelativeEntropyShannonDifferential{D <: DivergenceDefinition, E <: Entropy} <: Divergence
    e::E
    definition::D
    function RelativeEntropyShannonDifferential(; base = 2,
            definition::D = ShannonDivergenceDifferential()) where {D}
            e = Shannon(; base)
        new{D, typeof(e)}(e, definition)
    end
end
