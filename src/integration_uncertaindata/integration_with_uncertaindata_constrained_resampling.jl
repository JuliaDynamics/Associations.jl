
################################################################
# Integration with UncertainData.jl (with constraints)
################################################################

include("validate_constraints.jl")

"""
    causality(source, target, test::CausalityTest, resampling::ConstrainedResampling)

## Examples 

```julia
# Vectors of uncertain values. Also, convert to uncertain datasets
uvals_x = [UncertainValue(Normal, rand(Normal(0, 5)), abs(rand(Normal(0, 3)))) for i = 1:100]
uvals_y = [UncertainValue(Normal, rand(Normal(0, 5)), abs(rand(Normal(0, 3)))) for i = 1:100];
uvx = UncertainValueDataset(uvals_x)
uvy = UncertainValueDataset(uvals_y)

# Draw a single realisation of `uvx` and a single realisation of `uvy`
x, y = resample(uvx), resample(uvy)

# A transfer entropy test using the [`TransferOperatorGrid`](@ref) estimator.
test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = 1)

# `causality` accepts all combinations of real-valued and uncertain data types. Create some 
# constraints to apply when testing 
onevar_constraints = ConstrainedResampling(TruncateStd(1))
twovar_constraints = ConstrainedResampling(TruncateStd(2), TruncateStd(1))

# Perform causality tests
# ------------------------
# No constraints are needed, because no data are uncertain
causality(x, y, test_tog) 

# Only one variable is constrained, so provide only one set of constraints
causality(x, uvy, onevar_constraints, test_tog) 
causality(uvx, y, onevar_constraints, test_tog)
causality(x, uvals_y, onevar_constraints, test_tog)
causality(uvals_x, y, onevar_constraints, test_tog)

# Both variables are constrained, so provide two sets of constraints
causality(uvx, uvy, twovar_constraints, test_tog)
causality(uvals_x, uvals_y, twovar_constraints, test_tog) 
```
"""
function causality(source, target,  resampling::ConstrainedResampling, test::CausalityTest) end 

function causality(source, target, resampling::ConstrainedResampling{N}, test::TT) where {N, TT <: CausalityTest}

    validate_constraints(source, target, resampling)
    
    causality(source, target, test)
end


function causality(source::UVT1, target::UVT2, resampling::ConstrainedResampling{N}, test::TT) where {
            UVT1<:Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}},
            UVT2<:Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}},
            N,
            TT <: CausalityTest}
    
    validate_constraints(source, target, resampling)
    
    causality(resample(source, resampling[1]), resample(target, resampling[2]), test)
end

function causality(source, target::UVT, resampling::ConstrainedResampling{N}, test::TT) where {
            UVT <: Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}}, 
            N,
            TT <: CausalityTest}
    
    validate_constraints(source, target, resampling)

    causality(source, resample(target, resampling[1]), test)
end

function causality(source::UVT, target, resampling::ConstrainedResampling{N}, test::TT) where {
            UVT <: Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}}, 
            N,
            TT <: CausalityTest}

    validate_constraints(source, target, resampling)

    causality(resample(source, resampling[1]), target, test)
end