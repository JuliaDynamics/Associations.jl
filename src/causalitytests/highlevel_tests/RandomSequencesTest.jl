import UncertainData: RandomSequences

"""
    RandomSequencesTest(test::CausalityTest, sequences_resampling::RandomSequences)

A causality test where the `test` is applied to multiple independent draws of 
the datasets, where each draw from a randomly selected consecutive chunk of 
points.

## Examples 

Say we want do do a cross mapping test on a set of time series with `N = 200`
points. To get some ensemble statistics, we want to do an ensemble analysis 
on 150 randomly selected chunks of the time series with lengths ranging from 
`N-30` to `N-1`.

```julia
# A cross mapping test with default parameters (don't use default parameters in a real analysis!)
cm_test = CrossMappingTest()

# Combine with the random sequence subsamling
rtest = RandomSequencesTest(, RandomSequences(150, N-30:N-1))
```
"""
struct RandomSequencesTest{CT <: CausalityTest, RS <: RandomSequences}
    test::CT
    sequences_resampling::RS
end

export RandomSequencesTest