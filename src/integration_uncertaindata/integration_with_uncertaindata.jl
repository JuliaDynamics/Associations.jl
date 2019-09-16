
function causality(source::Vector{<:UVT1}, target::Vector{<:UVT2}, test::CT) where {
        UVT1<:AbstractUncertainValue, 
        UVT2<:AbstractUncertainValue,
        CT<:CausalityTest}
    causality(resample.(source), resample.(target), test)
end

function causality(source, target::Vector{<:UVT}, test::CT) where { 
        UVT<:AbstractUncertainValue,
        CT<:CausalityTest}
    causality(source, resample.(target), test)
end

function causality(source::Vector{<:UVT}, target, test::CT) where {
        UVT<:AbstractUncertainValue,
        CT<:CausalityTest}
    causality(resample.(source), target, test)
end

function causality(source::UVT1, target::UVT2, test::CT) where {
        UVT1<:AbstractUncertainValueDataset, 
        UVT2<:AbstractUncertainValueDataset,
        CT<:CausalityTest}
    causality(resample.(source), resample.(target), test)
end

function causality(source, target::UVT, test::CT) where { 
        UVT<:AbstractUncertainValueDataset,
        CT<:CausalityTest}
    causality(source, resample.(target), test)
end

function causality(source::UVT, target, test::CT) where {
        UVT<:AbstractUncertainValueDataset,
        CT<:CausalityTest}
    causality(resample.(source), target, test)
end