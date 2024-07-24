export PCMCI

"""
    PCMCI <: GraphAlgorithm
    PCMCI(
        test_unconditional = SurrogateAssociationTest(PearsonCorrelation(); nshuffles = 30),
        test_conditional= SurrogateAssociationTest(PartialCorrelation(); nshuffles = 30),
        τmax = 10,
        pmax = 10,
        qmax = 1,
    )

The PCMCI algorithm, using `test_unconditional` for pairwise tests and
`test_conditional` for conditional tests.

The maximum dimension of the condition is given by `pmax`, and `τmax` gives the maximum
embedding lag.
"""
Base.@kwdef struct PCMCI{U, C, L, P, Q} <: GraphAlgorithm
    test_unconditional::U = SurrogateAssociationTest(PearsonCorrelation(); nshuffles = 30)
    test_conditional::C = SurrogateAssociationTest(PartialCorrelation(); nshuffles = 30)
    τmax::L = 10
    pmax::P = 10
    qmax::Q = 2
end

include("parent_selection.jl")
