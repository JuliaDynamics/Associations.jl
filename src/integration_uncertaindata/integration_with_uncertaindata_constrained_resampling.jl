
################################################################
# Integration with UncertainData.jl (with constraints)
################################################################

include("validate_constraints.jl")

"""
    causality(source, target, test::CausalityTest, resampling::ConstrainedResampling{N})

Test for a causal influence from `source` to `target` using the provided causality `test`.
Both `source` and `target` might be uncertain data, which is resampled according to the 
`resampling` scheme.

This method uses the machinery from 
[`UncertainData.jl`](https://github.com/kahaaga/UncertainData.jl) to constrain the 
furnishing distributions of the uncertain data during resampling, which offers complete 
control over the resampling procedure.

## Arguments 

- **`source`**: The source. May be a real-valued vector, a `Vector{AbstractUncertainValue}`,
    or a `AbstractUncertainValueDataset`.

- **`target`**: The target. May be a real-valued vector, a `Vector{AbstractUncertainValue}`,
    or a `AbstractUncertainValueDataset`.

- **`resampling`**: A [`ConstrainedResampling`](@ref) instance. For one time series, 
    provide the `ConstrainedResampling` constructor with either a single sampling constraint 
    or a vector of sampling constraints. For two time series, provide the 
    `ConstrainedResampling` constructor with a tuple, where each element in the tuple 
    is either a single sampling constraint or a vector of sampling constraints. The 
    sampling constraints comprising the first element of the tuple is then mapped to the 
    uncertain values in `source`, and the second tuple is mapped to the elements of `target`.

- **`test`**: An instance of a [causality test](@ref causality_tests).

## Sampling constraints

See [`UncertainData.jl`](https://github.com/kahaaga/UncertainData.jl) for a list of possible 
constraints and how to construct them. The examples below show a few alternatives.
See also [`ConstrainedResampling`](@ref).

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

##############################
# Perform causality tests
##############################
# No constraints are needed, because no data are uncertain
causality(x, y, test_tog) 

# Applying the same constraint to all elements of the uncertain vectors
# ---------------------------------------------------------------------
onevar_constraints = ConstrainedResampling(TruncateStd(1))
twovar_constraints = ConstrainedResampling(TruncateStd(2), TruncateStd(1))


# Only one variable is constrained, so provide only one set of constraints
causality(x, uvy, test_tog, onevar_constraints) 
causality(uvx, y, test_tog, onevar_constraints)
causality(x, uvals_y, test_tog, onevar_constraints)
causality(uvals_x, y, test_tog, onevar_constraints)

# Both variables are constrained, so provide two sets of constraints
causality(uvx, uvy, test_tog, twovar_constraints)
causality(uvals_x, uvals_y, test_tog, twovar_constraints) 

# Applying different constraints on each element of the uncertain vectors
# -----------------------------------------------------------------------
constraints_x = [TruncateStd(1 + rand(Uniform(-0.5, 0.5))) for uv in uvals_x]
constraints_y = [TruncateStd(1 + rand(Uniform(-0.2, 0.2))) for uv in uvals_y]

# Constraints for x alone, constraints for y alone
crx = ConstrainedResampling(constraints_x)  
cry = ConstrainedResampling(constraints_y) 

# Compute causality statistic between one real-valued time series and an 
# uncertain time series. Resample the uncertain values according to the 
# `ConstrainedResampling`s we just defined.
causality(x, uvy, test_tog, cry) 
causality(uvx, y, test_tog, crx)
causality(x, uvals_y, test_tog, cry)
causality(uvals_x, y, test_tog, crx)

# Compute causality statistic between two uncertain time series. 
# We define separately constraints for y and y
crxy = ConstrainedResampling(constraints_x, constraints_y) 
causality(uvx, uvy, test_tog, crxy)
causality(uvals_x, uvals_y, test_tog, crxy) 
```
"""
function causality(source, target, test, resampling::ConstrainedResampling) end 

function causality(source, target, test::TT, resampling::ConstrainedResampling{N}) where {N, TT <: CausalityTest}

    validate_constraints(source, target, resampling)
    
    causality(source, target, test)
end


function causality(source::UVT1, target::UVT2, test::TT, resampling::ConstrainedResampling{N}) where {
            UVT1<:Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}},
            UVT2<:Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}},
            N,
            TT <: CausalityTest}
    
    validate_constraints(source, target, resampling)
    
    causality(resample(source, resampling[1]), resample(target, resampling[2]), test)
end

function causality(source, target::UVT, test::TT, resampling::ConstrainedResampling{N}) where {
            UVT <: Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}}, 
            N,
            TT <: CausalityTest}
    
    validate_constraints(source, target, resampling)

    causality(source, resample(target, resampling[1]), test)
end

function causality(source::UVT, target, test::TT, resampling::ConstrainedResampling{N}) where {
            UVT <: Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}}, 
            N,
            TT <: CausalityTest}

    validate_constraints(source, target, resampling)

    causality(resample(source, resampling[1]), target, test)
end