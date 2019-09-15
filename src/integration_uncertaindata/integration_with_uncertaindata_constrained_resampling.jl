
################################################################
# Integration with UncertainData.jl (with constraints)
################################################################

include("validate_constraints.jl")

function causality(source, target, test::TT, resampling::ConstrainedResampling{N}) where {N, TT <: CausalityTest}

    validate_constraints(source, target, resampling)
    
    causality(source, target, test)
end


function causality(source::UVT1, target::UVT2, test::TT, resampling::ConstrainedResampling{N}) where {
            UVT1<:Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}},
            UVT2<:Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}},
            TT <: CausalityTest,
            N}
    
    validate_constraints(source, target, resampling)
    
    causality(resample(source, resampling[1]), resample(target, resampling[2]), test)
end

function causality(source, target::UVT, test::TT, resampling::ConstrainedResampling{N}) where {
            UVT <: Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}}, 
            TT <: CausalityTest,
            N}
    
    validate_constraints(source, target, resampling)

    causality(source, resample(target, resampling[1]), test)
end

function causality(source::UVT, target, test::TT, resampling::ConstrainedResampling{N}) where {
            UVT <: Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}}, 
            TT <: CausalityTest,
            N}

    validate_constraints(source, target, resampling)

    causality(resample(source, resampling[1]), target, test)
end