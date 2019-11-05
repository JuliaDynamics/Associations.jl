import UncertainData.ConstrainedIndexValueResampling

"""
    ConstrainedTest(test::CausalityTest, constraints::ConstrainedResampling, n::Int)
    ConstrainedTest(test::CausalityTest, constraints::ConstrainedResampling)

A causality test where the supports of the uncertain values furnishing the 
uncertain values in the datasets are truncated (according to the provided 
`constraints`) prior to performing the `test`. 

`n` controls the number of independent resamplings over which the test 
is performed (if not provided, `n` is set to `1` by default).

## Examples

Assume we want to apply a causality test to two uncertain datasets `X` and `Y`,
but restricting the ranges of values their elements can take while resampling. 

```julia
# Constraints for X.indices and X.values
cx = (TruncateQuantiles(0.3, 0.7), TruncateStd(1.5))

# Constraints for Y.indices and Y.values
cy = (TruncateQuantiles(0.3, 0.7), TruncateStd(1.5))

# Need to gather to feed to the ConstrainedTest constructor.
cs = ConstrainedIndexValueResampling(cx, cy)

# A cross mapping test applied to 100 independent realisations
# of `X` and `Y`, resampled after constraining the data.
ctest = ConstrainedTest(CrossMappingTest(), cs, 100)
```
"""
struct ConstrainedTest{
        CT <: CausalityTest, 
        CR <: ConstrainedIndexValueResampling{N_VARIABLES, N_DATASETS} 
            where {N_VARIABLES, N_DATASETS}} <: PreprocessedDataCausalityTest
    test::CT
    constrained_resampling::CR
    n::Int
end

function ConstrainedTest(test, constrained_resampling)
    ConstrainedTest(test, constrained_resampling, 1)
end

export ConstrainedTest