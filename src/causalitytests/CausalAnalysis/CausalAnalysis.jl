

""" 
    CausalAnalysis(test, data, result)

The result of a causal analysis.

## Fields

- **`test`**: The parameters of the test.
- **`data`**: The data on which the test was applied.
- **`result`**: The result of applying the relevant causality statistic using
the parameters in `test` to `data`. 

## Implements 

- [`summarise(f::Function, analysis::CausalAnalysis, args...; kwargs...)`]. Summarise 
a causal analysis using some summary statistic.
"""
struct CausalAnalysis{CT <: Union{CausalityTest, MetaCausalityTest}, DT, RT} <: AbstractCausalAnalysis{CT, DT, RT}
    test::CT
    data::DT
    result::RT
end

function Base.show(io::IO, x::AbstractCausalAnalysis)
    T = typeof(x).name

    CT = typeof(x.test)
    RT = typeof(x.result)
    DT = typeof(x.data)
    CT_str = "\n  test: $CT"
    DT_str = "\n  data: $DT"
    RT_str = "\n  result: $(RT)"
    summary = join([CT_str, DT_str, RT_str], "")
    print(io, "$(T)$summary\n")
end

"""
    causality(test::Union{CausalityTest, MetaCausalityTest}, x, y)
    causality(test::Union{CausalityTest, MetaCausalityTest}, x, y, z)

Apply the causality `test` to the provided data series. Providing the 
test as the first argument returns a [`CausalAnalysis`](@ref) 
instance which stores the test, the data used for the test, and the 
result in a single struct.
"""
function causality(test::Union{CausalityTest, MetaCausalityTest}) end

function causality(test::Union{CausalityTest, MetaCausalityTest}, x, y)
    CausalAnalysis(test, (x, y), causality(x, y, test))
end

function causality(test::Union{CausalityTest, MetaCausalityTest}, x, y, z)
    CausalAnalysis(test, (x, y, z), causality(x, y, z, test))
end
