function validate_constraints(
        source::UVT1, 
        target::UVT2, 
        resampling::ConstrainedResampling{N}) where {
            UVT1<:Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}}, 
            UVT2<:Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}},
            N}
    if N != 2
        msg = "Both `source` and `target` are uncertain, so `resampling` must contain exactly two separate sets of sampling constraints"    
        throw(ArgumentError(msg))
    end
end

function validate_constraints(
        source::Vector{<:Real}, 
        target::UVT, 
        resampling::ConstrainedResampling{N}) where {
            UVT<:Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}}, N}

    if N != 1
        msg = "Only `target` is uncertain, but $N sets of constraints were provided. Only one should be given."
        throw(ArgumentError(msg))
    end
end

function validate_constraints(
        source::UVT, 
        target::Vector{<:Real}, 
        resampling::ConstrainedResampling{N}) where {
            UVT<:Union{AbstractUncertainValueDataset, Vector{<:AbstractUncertainValue}}, N}

    if N != 1
        msg = "Only `source` is uncertain, but $N sets of constraints were provided. Only one should be given."
        throw(ArgumentError(msg))
    end
end

function validate_constraints(
        source::Vector{<:Real}, 
        target::Vector{<:Real}, 
        resampling::ConstrainedResampling{N}) where {
            UVT<:AbstractUncertainValue, N}
    @warn "Neither `source` nor `target` are uncertain. Not constraining them."
end
