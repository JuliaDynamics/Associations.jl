export ShannonDivergence

"""
    ShannonDivergence  <: DivergenceDefinition
    ShannonDivergence()

An instruction to compute the Shannon relative entropy. Also called the KL divergence.

[^Bulinski2021]:
    Bulinski, A., & Dimitrov, D. (2021). Statistical estimation of the Kullback-Leibler
    divergence. Mathematics, 9(5), 544.
"""
struct ShannonDivergence <: DivergenceDefinition end
