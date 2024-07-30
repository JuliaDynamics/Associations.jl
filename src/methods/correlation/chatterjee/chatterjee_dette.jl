using Random

export ChatterjeeDetteCorrelation
export Chatterjee


# Based on the description in https://arxiv.org/pdf/2008.11619
"""
    ChatterjeeDetteCorrelation <: CorrelationMeasure
    ChatterjeeDetteCorrelation(x::AbstractVector, y::AbstractVector; 
        jitter = false, handle_ties = false,
        rng = Random.default_rng(), 
        lt = (a, b) -> isless_rand(rng, a, b))

The Chatterjee-Dette correlation measure is an asymmetric measure of dependence between
two variables.

The input data `x` and `y` are required, because they are used to pre-allocate container for 
handle_tieser computations.

## Description 

The correlation measure is defined as

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

!!! info "Speeding up computations"
    If `handle_ties == true`, then the first formula is used. Use this 
    if you know for sure that the input data has no ties. If 
    `handle_ties == false`, then the second formula is used (
    which is faster).
    
    For numerical input data, `jitter == true` will set `handle_ties == true` and add a 
    small amount of uniform noise to the data. This will almost almost surely break 
    any ties before the second formula is applied.
"""
struct ChatterjeeDetteCorrelation{
        I <: Integer, 
        J <: Integer, 
        V <: Union{AbstractVector, Nothing}
    } <: CorrelationMeasure

    # A copy of the original data, to which we add noise.
    xc::V

    # A container holding the indices that would sort the *first* input vector (X).
    x_sortidxs::Vector{I}

    # The ranks. r[i] := sum_{j=1}^n I(Y_j ≤ Y_i), where Y is the *second* input variable.
    rs::Vector{J}

    jitter::Bool

    handle_ties::Bool


    function ChatterjeeDetteCorrelation(x, y; jitter = false, handle_ties = true)
        @assert length(x) == length(y)
        n = length(x)
        x_sortidxs = zeros(Int, n)
        rs = zeros(Int, n)

        # If we jitter, we can almost surely guarantee that there are no ties and use
        # the handle_ties version.
        if jitter  && x isa AbstractVector{<:Number}
            if x isa AbstractVector{<:Number}
                xc = zeros(Float64, n)
                new{Int, Int, Nothing}(xc, x_sortidxs, rs, jitter)
            else
                msg = "To add jitter to the input data, the input data must be a numeric vector"
                throw(ArgumentError(msg))
            end
        else
            new{Int, Int, Nothing}(nothing, x_sortidxs, rs, jitter)
        end
    end
end

struct Chatterjee{XS, R, L, DR, YR, YSI, YSS, RNG} <: CorrelationMeasure
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
function Chatterjee(; handle_ties = true, rng = Random.default_rng(), lt = (a, b) -> isless_rand(rng, a, b))
    return Chatterjee{Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, typeof(rng)}(
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, rng, handle_ties, lt)
end

# Preallocated version
function Chatterjee(x, y; handle_ties = true, rng = Random.default_rng(), lt = (a, b) -> isless_rand(rng, a, b))
    verify_inputs_chatterjee(x, y)
    n = length(x)
    xs_inds = zeros(Int, n)
    y_sortedbyx = similar(y)  # rearranged `y` (sorted according to `xs_inds`)
    y_sortedbyx_inds = zeros(Int, n)  # rearranged `y` (sorted according to `xs_inds`)
    y_sortedbyx_sorted = similar(y) # sort the rearranged `y`.

    𝓁 = zeros(Int, n)
    r = zeros(Int, n)
    Δr = zeros(n)
    return Chatterjee(xs_inds, r, 𝓁, Δr, y_sortedbyx, y_sortedbyx_inds, y_sortedbyx_sorted, rng, handle_ties, lt)
end

function association(m::Chatterjee{Nothing}, x, y)
    verify_inputs_chatterjee(x, y)
    n = length(x)

    y_sortedbyx = y[sortperm(x, lt = m.lt)] # Rearranged `y`s
    if m.handle_ties
        r, 𝓁 = loopcount_rs_ℓs(y_sortedbyx)
    else
        r = loopcount_rs(y_sortedbyx)
        𝓁 = nothing
    end
    
    # Absolute values of the first differences for the ranks.
    Δr = zeros(Int, n - 1)
    Δr = absfirstdiffs!(Δr, r)

    chatterjee_coefficient(Δr, 𝓁, n, m.handle_ties)
end

# Pre-allocated version
function association(m::Chatterjee{<:AbstractVector}, x, y)
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

function loopcount_rs(y_sortedbyx)
    sorted_y_sortedbyx = sort(y_sortedbyx)
    r = zeros(length(sorted_y_sortedbyx))
    loopcount_rs!(r, sorted_y_sortedbyx, y_sortedbyx)
    return r
end

function loopcount_rs!(r, sorted_y_sortedbyx, y_sortedbyx)    
    for (i, yᵢ) in enumerate(y_sortedbyx)
        # The position of the last (rightmost) occurrence of yᵢ in the 
        # sorted data is the rank. 
        rightmost_occurrence = searchsortedlast(sorted_y_sortedbyx, yᵢ)
        r[i] = rightmost_occurrence
    end
end

function loopcount_rs_ℓs(y_sortedbyx)
    n = length(y_sortedbyx)
    sorted_y_sortedbyx = sort(y_sortedbyx)
    r = zeros(Int, n)
    𝓁 = zeros(Int, n)
    loopcount_rs_ℓs!(r, 𝓁, sorted_y_sortedbyx, y_sortedbyx)
    return r, 𝓁
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

Compute the Chatterjee-Dette correlation coefficient given the length of the 
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
        msg = "`length(x)` must equal `length(y)` for `ChatterjeeDetteCorrelation`. " * 
            "Got length(x)=$(length(x)) and length(y)=$(length(x))".
        throw(ArgumentError(msg))
    end
end

function isless_rand(rng, a, b)
    if  a < b
        true
    elseif a > b
        false
    else
        rand(rng, Bool)
    end
end

