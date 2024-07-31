using Random

export ChatterjeeCorrelation


# Based on the description in https://arxiv.org/pdf/2008.11619
"""
    ChatterjeeCorrelation <: CorrelationMeasure
    ChatterjeeCorrelation(; handle_ties = true, rng = Random.default_rng())

The Chatterjee correlation measure is an asymmetric measure of dependence between
two variables. 

!!! info "Speeding up computations"
    If `handle_ties == true`, then the first formula below is used. If you know
    for sure that there are no ties in your data, then set `handle_ties == false`, 
    which will use the second (faster) formula below.

!!! note "Randomization and reproducibility"
    When rearranging the input datasets, the second variable `y` is sorted 
    according to a sorting of the first variable `x`. If `x` has ties, then 
    these ties are broken randomly and uniformly. For complete reproducibility in 
    this step, you can specify `rng`. If `x` has no ties, then no randomization 
    is performed.

# Usage

- Use with [`association`](@ref) to compute the raw Chatterjee correlation coefficient.
- Use with [`SurrogateAssociationTest`](@ref) to perform a surrogate test for significance 
    of a Chatterjee-type association ([example](@ref example_SurrogateAssociationTest_ChatterjeeCorrelation)). 
    When using a surrogate test for significance, the *first* input variable is shuffled 
    according to the given surrogate method.

## Description 

The correlation statistic is defined as

```math
\\epsilon_n(X, Y) = 
1 - \\dfrac{n\\sum_{i=1}^{n-1} |r_{i+1} - r_i|}{2\\sum_{i=1}^n }.
```

When there are no ties among the ``Y_1, Y_2, \\ldots, Y_n``, the 
measure is 

```math
\\epsilon_n(X, Y) = 
1 - \\dfrac{3\\sum_{i=1}^{n-1} |r_{i+1} - r_i|}{n^2 - 1}.
```

This statistic estimates a quantity proposed by [Dette2013](@citet), as indicated in 
[Shi2022](@citet). It can therefore also be called the Chatterjee-Dette-Siburg-Stoimenov
correlation coefficient.

## Estimation

- [Example 1](@ref example_ChatterjeeCorrelation). Estimating the Chatterjee correlation coefficient
    for independent and for dependent variables. 
- [Example 2](@ref example_SurrogateAssociationTest_ChatterjeeCorrelation). Testing the significance
    of a Chatterjee-type association using a surrogate test.
"""
struct ChatterjeeCorrelation{XS, R, L, DR, YR, YSI, YSS, RNG} <: CorrelationMeasure
    # Pre-allocated containers (if input data are given).
    # If no input data are given, then these are all `nothing`.
    xs_inds::XS # indices that would sort `x`
    r::R # ranks
    𝓁::L # 𝓁s
    Δr::DR # absolute value of first differences of ranks
    y_sortedbyx::YR # ys, sorted according to the order that sorts the xs.
    y_sortedbyx_inds::YSI # A container holding the indices that would sort the *second* input vector (Y).
    y_sortedbyx_sorted::YSS # the sorted version of `y_sortedbyx`
    rng::RNG
    handle_ties::Bool
    lt::Function
end

# Non-preallocated version
function ChatterjeeCorrelation(; handle_ties = true, rng = Random.default_rng(), lt = (a, b) -> isless_rand(rng, a, b))
    return ChatterjeeCorrelation{Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, typeof(rng)}(
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, rng, handle_ties, lt)
end

# Preallocated version
function ChatterjeeCorrelation(x, y; 
        handle_ties = true, rng = Random.default_rng(), lt = (a, b) -> isless_rand(rng, a, b))
    verify_inputs_chatterjee(x, y)
    n = length(x)
    xs_inds = zeros(Int, n)
    y_sortedbyx = similar(y)  # rearranged `y` (sorted according to `xs_inds`)
    y_sortedbyx_inds = zeros(Int, n)  # rearranged `y` (sorted according to `xs_inds`)
    y_sortedbyx_sorted = similar(y) # sort the rearranged `y`.

    𝓁 = zeros(Int, n)
    r = zeros(Int, n)
    Δr = zeros(n)
    return ChatterjeeCorrelation(xs_inds, r, 𝓁, Δr, y_sortedbyx, y_sortedbyx_inds, 
        y_sortedbyx_sorted, rng, handle_ties, lt)
end

# If the user provided a non-preallocated version, then we simply create a new version
# of the measure that *is* pre-allocated and use that.
function association(m::ChatterjeeCorrelation{Nothing}, x, y)
    measure = ChatterjeeCorrelation(x, y, rng = m.rng, lt = m.lt, handle_ties = m.handle_ties)
    return association(measure, x, y)
end

# Pre-allocated version
function association(m::ChatterjeeCorrelation{<:AbstractVector}, x, y)
    verify_inputs_chatterjee(x, y)
    n = length(x)

    # Rearrange input data according to paper.
    m.xs_inds[:] = sortperm!(m.xs_inds, x, lt = m.lt) 
    m.y_sortedbyx[:] = y[m.xs_inds]

    # Sort the rearranged `y`s so that we can more efficiently find the ranks.
    m.y_sortedbyx_inds[:] = sortperm!(m.y_sortedbyx_inds, m.y_sortedbyx)
    m.y_sortedbyx_sorted[:] = m.y_sortedbyx[m.y_sortedbyx_inds]
    
    # If we explicitly deal with ties, then we need to find both the ranks `r`s and 
    # the `𝓁`s that go into the denominator of the final expression. If we don't deal 
    # with ties, then we only need to find the ranks `r`, which save some time.
    if m.handle_ties
        loopcount_rs_ℓs!(m.r, m.𝓁, m.y_sortedbyx_sorted, m.y_sortedbyx);
    else
        loopcount_rs!(m.r, m.y_sortedbyx_sorted, m.y_sortedbyx)
    end
    
    # First differences of the ranks, stored in  the pre-allocated `Δr`.`
    absfirstdiffs!(m.Δr, m.r)

    return chatterjee_coefficient(m.Δr, m.𝓁, n, m.handle_ties)
end

function loopcount_rs!(r, sorted_y_sortedbyx, y_sortedbyx)    
    for (i, yᵢ) in enumerate(y_sortedbyx)
        # The position of the last (rightmost) occurrence of yᵢ in the 
        # sorted data is the rank. 
        rightmost_occurrence = searchsortedlast(sorted_y_sortedbyx, yᵢ)
        r[i] = rightmost_occurrence
    end
end

function loopcount_rs_ℓs!(r, ℓ, sorted_y_sortedbyx, y_sortedbyx)
    n = length(y_sortedbyx)
    
    for (i, yᵢ) in enumerate(y_sortedbyx)
        # Indices of the first and last occurrences of `yᵢ` in the sorted data.
        leftmost_occurrence::Int = searchsortedfirst(sorted_y_sortedbyx, yᵢ)
        rightmost_occurrence::Int = searchsortedlast(sorted_y_sortedbyx, yᵢ)
        
        # rᵢ (rank of yᵢ): count of elements less than or equal to yᵢ
        r[i] = rightmost_occurrence
        
        # ℓ: count of elements greater than or equal to yᵢ
        ℓ[i] = n - leftmost_occurrence + 1
    end
    
    return r, ℓ
end

"""
    chatterjee_coefficient(Δr, 𝓁, n::Integer, handle_ties::Bool)

Compute the Chatterjee correlation coefficient given the length of the 
input `n`, pre-computed first differences of the ranks `Δr` and pre-computed
`𝓁`.
"""
function chatterjee_coefficient(Δr, 𝓁, n::Integer, handle_ties::Bool)
    if handle_ties
        num = n*sum(Δr)
        den = 2*sum(𝓁ᵢ * (n - 𝓁ᵢ) for 𝓁ᵢ in 𝓁)
    else
        num = 3*sum(Δr)
        den = (n^2 - 1)
    end
    return 1 - num / den
end

"""
    absfirstdiffs!(dx, x)

Compute absolute values of the first differences between the elements of `x` into the 
preallocated x::AbstractVectoray `dx` where `length(dx) == length(x)` (not checked).
"""
function absfirstdiffs!(dx, x)
    n = length(x)
    for i in 1:(n - 1)
        dx[i] = abs(x[i + 1] - x[i])
    end
    return dx
end

function verify_inputs_chatterjee(x, y)
    if length(x) != length(y)
        msg = "`length(x)` must equal `length(y)` for `ChatterjeeCorrelation`. " * 
            "Got length(x)=$(length(x)) and length(y)=$(length(x))".
        throw(ArgumentError(msg))
    end
end

# break ties randomly.
function isless_rand(rng, a, b)
    if  a < b
        true
    elseif a > b
        false
    else
        rand(rng, Bool)
    end
end

