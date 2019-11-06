import UncertainData:
    AbstractUncertainValue,
    AbstractUncertainValueDataset, 
    AbstractUncertainIndexValueDataset,
    resample

"""
    resample_and_subset(x, r::AbstractVector{Int})

Resample `x` assuming the points in `x` are independent,
then return the points with indices `r`.
"""
function resample_and_subset(x, r::AbstractVector{Int})
    v = zeros(Float64, length(r))
    
    if x isa AbstractVector
        v[:] = x[r]
    elseif typeof(x) <: Union{Vector{<:AbstractUncertainValue}, AbstractUncertainValueDataset}
        v[:] = resample(x)[r]
    elseif typeof(x) <: AbstractUncertainIndexValueDataset
        v[:] = resample(x)[2][r]
    end
    
    return v
end

"""
    causality(x, y, test::RandomSequencesTest)

Apply a causality test on random, consecutive chunks of `x` and `y`.

The chunk length can be fixed (an integer) or a collection of
chunk lengths.

Works on any inputs `x` and `y`, both scalar vectors and 
uncertain datasets.

## Example 

For scalar valued vectors:

```julia 
x, y = rand(100), rand(100)

# Use a cross-mapping test with default parameters (not recommended,
# make sure you understand the parameters!)
cm_test = CrossMappingTest()

# A random sequences test that applies the cross mapping test 
# to 150 different random chunks of the data of lengths 80 to 90
n_realisations, chunk_lengths = 150, 80:90
rs = RandomSequences(n_realisations, chunk_lengths)

# Combine to a RandomSequencesTest
rs_test = RandomSequencesTest(cm_test, chunk_lengths)

# Compute causality statistic
causality(x, y, rs_test)
```

This also works on uncertain data, or any combination of scalar vectors
and uncertain data:

```julia
# Some example data
N = 300
sys = ar1_unidir(c_xy = 0.8)
X, Y = example_uncertain_indexvalue_datasets(sys, N, (1, 2),
    d_xval = Uniform(0.001, 0.05), d_yval = Uniform(0.001, 0.05));

# Apply cross mapping test on 100 different randomly selected chunks
# of lengths 75 to 90
r = RandomSequences(100, 75:90)
rs_test = RandomSequencesTest(CrossMappingTest(), r)
causality(X, Y, rs_test)
```
"""
function causality(
        x::Union{AbstractVector, Vector{<:AbstractUncertainValue}, AbstractUncertainValueDataset, AbstractUncertainIndexValueDataset}, 
        y::Union{AbstractVector, Vector{<:AbstractUncertainValue}, AbstractUncertainValueDataset, AbstractUncertainIndexValueDataset},
        test::RandomSequencesTest)
    
    seq_length = test.sequences_resampling.sequence_length
    n = test.sequences_resampling.n
    N = length(x)
    results = [Vector{Float64}(undef, 0) for i = 1:n]
    
    for i = 1:n
        if seq_length isa Int
            idx_start = rand(1:(length(x) - seq_length))
            r = idx_start:(idx_start + seq_length)
        elseif seq_length isa AbstractVector{Int}
            # seq_length is now something we can sample directly from
            seqlen = rand(seq_length)
            idx_start = rand(1:(N + 1 - seqlen))
            r = idx_start:(idx_start + seqlen - 1)
        else
            throw(ArgumentError("`resampling.sequence_length`must be an integer or a collection of integers"))
        end
        xs = resample_and_subset(x, r)
        ys = resample_and_subset(y, r)
        res = causality(xs, ys, test.test)
        
        append!(results[i], res)
    end
    
    return results
end

export causality, resample_and_subset