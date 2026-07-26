
export PCMCI
export PCMCISelectedParents
export PCMCIParent

"""
    PCMCI <: GraphAlgorithm
    PCMCI(; τmax::Int = 5, pmax::Int = -1, qmax::Int = 1, pX::Int = 3, 
        αPC = 0.15, α = 0.05,
        fdr_adjust = true,
        utest = SurrogateAssociationTest(KSG2(MIShannon(), k = 10), nshuffles = 19),
        ctest = LocalPermutationTest(MesnerShalizi(CMIShannon(), k=10, w=0), nshuffles = 19),
    )

The PCMCI causal network inference algorithm [Runge2019](@cite). 

# Keyword arguments

The PCMCI algorithm is here implemented with compatibility with *any* of the pairwise and 
conditional independence tests implemented in Associations.jl.

- **`utest`**: The pairwise/unconditional association test. Defaults to a [`SurrogateAssociationTest`](@ref)
    with a mutual information estimator, but can be any [`IndependenceTest`](@ref).
- **`ctest`**: The conditional association test. Defaults to a [`SurrogateAssociationTest`](@ref)
    with a conditional mutual information estimator, but can be any conditional [`IndependenceTest`](@ref).
- **`αPC`**: PC1 removal threshold (regularization parameter). Too small `αPC` increases false positives,
    while too small `αPC` reduces detection power and increases runtime.
- **`α`**: Significance threshold in the MCI step. 
- **`τmax`**: The maximum lag to consider. 
- **`pmax`**: The maximum dimension of the condition. Use `-1` (the default) for no
    restriction, corresponding to the paper's `pmax = Nτmax`, where `N` is the number of
    variables (`dimension(input_dataset)`) and `τmax` is the maximum lag.
- **`qmax`**: The maximum number of combinations.
- **`pX`**: The maximum number of the driver variable's strongest parents to condition on in the MCI stage 
    (last term of Eq. 3 in [Runge2019](@cite)). Conditioning on these accounts for autocorrelation. Too low 
    risks `pX` inflated false positives, too high lowers test power and increases dimensionality (and hence 
    increases computation time).
- **`fdr_adjust`**: Adjust p-values for all links using the False Discovery Rate (FDR) approach?
- **`use_abs_utest`**: Whether to rank parents by the absolute value `|I|` of the *pairwise*
    (`utest`) test statistic (`true`), or by the raw value (`false`). Defaults to `true`.
- **`use_abs_ctest`**: Whether to rank parents by the absolute value `|I|` of the *conditional*
    (`ctest`) test statistic (`true`), or by the raw value (`false`). Defaults to `true`.

    Using `|I|` (the paper's convention) is meaningful for signed measures like partial
    correlation, where `I = -0.8` is a *stronger* dependence than `I = +0.3`. For unsigned
    measures like MI/CMI, some estimators may yield negative values despite a theoretical
    minimum of `0`; in these cases, set the corresponding flag to `false` to compare the raw
    values directly. The flags are separate because `utest` and `ctest` may use different
    measures with different sign behaviour.

## Used with 

- [`infer_graph`](@ref). The input variables must be given as a `StateSpaceSet` or a 
    `Vector{Vector{<:Real}}`.
## Returns

When used with [`infer_graph`](@ref), it returns a `Vector{`[`PCMCISelectedParents`](@ref)`}`,
where `p[i]` are the parents for the variable `Xⁱₜ`. This result can be converted to a
`SimpleDiGraph` from Graphs.jl by calling `SimpleDiGraph(res)`.
"""
Base.@kwdef struct PCMCI{U,C,T} <: GraphAlgorithm
    utest::U = SurrogateAssociationTest(KSG2(MIShannon(), k=10, w=0), nshuffles=19)
    ctest::C = LocalPermutationTest(MesnerShalizi(CMIShannon(), k=10, w=0), nshuffles=19)
    τmax::T = 5
    pmax::Int = -1 # maximum condition dimension; -1 = unrestricted (N*τmax)
    qmax::Int = 1 # maximum number of combinations
    pX::Int = 3
    αPC = 0.2
    α = 0.05
    fdr_adjust = true
    use_abs_utest::Bool = true
    use_abs_ctest::Bool = true
end


"""
    PCMCIParent(i, τ, pval, test_statistic)

The index `i` and lag `τ` of a parent node, along with the p-value and 
test statistic value for the link.
"""
struct PCMCIParent
    i::Int
    τ::Int
    pval::Real
    test_statistic::Real
end
"""
    PCMCISelectedParents <: AbstractSelectedParents

Stores the parents of a single variable `xᵢ(0)` inferred by the [`PCMCI`](@ref) algorithm.
When [`PCMCI`](@ref) is used with [`infer_graph`](@ref), a `Vector{PCMCISelectedParents}`
is returned, i.e. one element per input variable, where `p[i]` holds the parents of `xᵢ(0)`.

## Fields

- `i`: index of the target variable (`xᵢ(0)` is the target).
- `parents_js`: variable index of each selected parent.
- `parents_τs`: lag of each selected parent.
- `pvals`: the (possibly FDR-adjusted) p-value of each selected link.
- `test_statistics`: the test-statistic value of each selected link.

When printed in the console the `[p=…, I=…]` values are displayed for each parent.
"""
struct PCMCISelectedParents{PJ,PT,PV,TS} <: AbstractSelectedParents
    i::Int
    parents_js::PJ
    parents_τs::PT
    pvals::PV
    test_statistics::TS
end

# Convenience constructor from a vector of per-link `PCMCIParent`s.
function PCMCISelectedParents(i::Int, links::AbstractVector{PCMCIParent})
    return PCMCISelectedParents(
        i,
        [l.i for l in links],
        [l.τ for l in links],
        [l.pval for l in links],
        [l.test_statistic for l in links],
    )
end

# Include the p-value and test statistic after each parent when printing (see
# `print_condvars(::AbstractSelectedParents)`).
function parent_annotation(p::PCMCISelectedParents, k::Int)
    return " [p=$(round(p.pvals[k]; digits=4)), I=$(round(p.test_statistics[k]; digits=4))]"
end
