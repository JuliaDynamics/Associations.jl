using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: ContinuousDynamicalSystem
using DynamicalSystemsBase: trajectory
using SimpleDiffEq: SimpleATsit5
using Distributions: Uniform

export RosslerForced9

Base.@kwdef struct RosslerForced9{V,W1,W2,W3,CXY,CYX,CZX,CZY,A1,A2,A3,B1,B2,B3,C1,C2,C3} <: ContinuousDefinition
    xi::V = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    ω₁::W1 = 1.015
    ω₂::W2 = 0.985
    ω₃::W3 = 0.95
    c_xy::CXY = 0.1
    c_yx::CYX = 0.1
    c_zx::CZX = 0.05
    c_zy::CZY = 0.05
    a₁::A1 = 0.15
    a₂::A2 = 0.2
    a₃::A3 = 10
    b₁::B1 = 0.15
    b₂::B2 = 0.2
    b₃::B3 = 10
    c₁::C1 = 0.15
    c₂::C2 = 0.2
    c₃::C3 = 10
end


function system(definition::RosslerForced9)
    return ContinuousDynamicalSystem(eom_rosslerforced9, definition.xi, definition)
end

@inline @inbounds function eom_rosslerforced9(u, p, t)
    (; xi, ω₁, ω₂, ω₃, c_xy, c_yx, c_zx, c_zy, a₁, a₂, a₃, b₁, b₂, b₃, a₃, c₁, c₂, c₃) = p
    x1, x2, x3, y1, y2, y3, z1, z2, z3 = u

    dx1 = -ω₁*(x2 + x3) + c_yx*(y1 - x1) + c_zx*(z1 - x1)
    dx2 = ω₁*x1 + a₁*x2
    dx3 = a₂ + x3*(x1 - a₃)

    dy1 = -ω₂*(y2 + y3) + c_xy*(x1 - y1) + c_zy*(z1 - y1)
    dy2 = ω₂*y1 + b₁*y2
    dy3 = b₂ + y3*(y1 - b₃)

    dz1 = -ω₂*(z2 + z3)
    dz2 = ω₂*z1 + c₁*z2
    dz3 = c₂ + z3*(z1 - c₃)

    return SVector{9}(dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3)
end
