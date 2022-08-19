export EffortToCompress, EffortToCompressSlidingWindow

"""
    EffortToCompress(; normalize = false, alphabet_size = nothing)

The effort-to-compress (ETC)[^Nagaraj2013] algorithm quantifies the compression complexity 
of a time series.

## Normalization

If `normalize == false`, then computes the number of compression steps `N` it takes for the 
symbol sequence to reach zero entropy (either a constant sequence, or a length-1 sequence).

For a length-`N` sequence, the maximum number of possible compression steps is `N-1`.
If `normalize == true`, then the ETC value is normalized to `N-1`, yielding a number on 
`[0, 1]`, where `0` indicates minimal compression complexity and `1` indicates maximal 
compression complexity. 

## Alphabet size

The `alphabet_size` is the number of possible values the time series can take. For example,
a binary time series, `alphabet_size = 2`. For a five-box binned time series, 
`alphabet_size = 5`. 

The `alphabet_size` parameter must be given a positive integer value in order to compute 
ETC for multivariate time series. For univariate time series, `alphabet_size` is ignored.

[^Nagaraj2013]: Nagaraj, N., Balasubramanian, K., & Dey, S. (2013). A new complexity 
    measure for time series analysis and classification. The European Physical Journal 
    Special Topics, 222(3), 847-860.
"""
struct EffortToCompress
    alphabet_size::Union{Nothing, Int}
    normalize::Bool

    function EffortToCompress(; 
        alphabet_size::Union{Nothing, Int} = nothing,
        normalize = false)
        return new(alphabet_size, normalize)
    end
end

"""
    EffortToCompressSlidingWindow(; normalize::bool = false, window_size = 15, step = 1)

Like [`EffortToCompress`](@ref), but applied on a constant-length sliding window across 
the time series, using the given `step` and `window_size`.
"""
struct EffortToCompressSlidingWindow
    alphabet_size::Union{Nothing, Int}
    normalize::Bool
    window_size::Int
    step::Int
    EffortToCompressSlidingWindow(; 
        alphabet_size::Union{Nothing, Int} = nothing,
        normalize::Bool = false,
        window_size::Int = 15,
        step::Int = 1,
        ) = new(alphabet_size, normalize, window_size, step)
end


include("etc_utils.jl")
include("etc_univariate.jl")
include("etc_multivariate.jl")
include("etc_joint.jl")
