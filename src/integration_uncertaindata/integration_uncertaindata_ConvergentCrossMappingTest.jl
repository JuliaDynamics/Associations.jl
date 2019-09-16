import ..CausalityTests: ConvergentCrossMappingTest

################################################################
# Integration with UncertainData.jl
################################################################
function causality(source::Vector{<:UVT1}, target::Vector{<:UVT2}, p::ConvergentCrossMappingTest) where {
        UVT1<:AbstractUncertainValue, 
        UVT2<:AbstractUncertainValue}
    causality(resample.(source), resample.(target), p)
end


function causality(source, target::Vector{<:UVT}, p::ConvergentCrossMappingTest) where { 
        UVT<:AbstractUncertainValue}
    causality(source, resample.(target), p)
end


function causality(source::Vector{<:UVT}, target, p::ConvergentCrossMappingTest) where {
        UVT<:AbstractUncertainValue}
    causality(resample.(source), target, p)
end


function causality(source::UVT1, target::UVT2, p::ConvergentCrossMappingTest) where {
        UVT1<:AbstractUncertainValueDataset, 
        UVT2<:AbstractUncertainValueDataset}
    causality(resample.(source), resample.(target), p)
end


function causality(source, target::UVT, p::ConvergentCrossMappingTest) where { 
        UVT<:AbstractUncertainValueDataset}
    causality(source, resample.(target), p)
end


function causality(source::UVT, target, p::ConvergentCrossMappingTest) where {
        UVT<:AbstractUncertainValueDataset}
    causality(resample.(source), target, p)
end

