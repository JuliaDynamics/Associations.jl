


"""
     AbstractCausalAnalysis{CT <: Union{CausalityTest, MetaCausalityTest}, DT, RT}

Abstract return type for a causality test. `CT` is the type of the causality test,
`DT` is the data type, and `RT` is the result type.

As a minimum, concrete subtypes must implement:

- **`test`**: The parameters of the test.
- **`data`**: The data on which the test was applied.
- **`result`**: The result of applying the relevant causality statistic using
    the parameters in `test` to `data`. 
"""
abstract type AbstractCausalAnalysis{CT <: Union{CausalityTest, MetaCausalityTest}, DT, RT} end
