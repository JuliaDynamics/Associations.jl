# This is a workaround until we can make faster horizontal concatenation in StateSpaceSets.jl
import StateSpaceSets: StateSpaceSet
using StaticArrays: SVector
export StateSpaceSet

function StateSpaceSet(x::VectorOrStateSpaceSet...; container = SVector)
    xs = (xᵢ isa AbstractStateSpaceSet ? Matrix(xᵢ) : reshape(xᵢ, length(xᵢ), 1) for xᵢ in x)
    StateSpaceSet(hcat(xs...); container)
end