using DynamicalSystemsBase
using Random
rng = Random.MersenneTwister(1234)

Base.@kwdef struct Logistic2Unidir{V, C, R1, R2, Σy, R}
    xi::V = [0.5, 0.5]
    c_xy::C = 0.1
    r₁::R1 = 3.78
    r₂::R2 = 3.66
    σ_xy::Σy = 0.05
    rng::R = Random.default_rng()
end

function eom_logistic2uni(u, p::Logistic2Unidir, t)
    (; xi, c_xy, r₁, r₂, σ_xy, rng) = p
    x, y = u
    f_xy = (y +  (c_xy*(x + σ_xy * rand(rng))/2) ) / (1 + (c_xy/2)*(1+σ_xy))

    dx = r₁ * x * (1 - x)
    dy = r₂ * (f_xy) * (1 - f_xy)
    return SVector{2}(dx, dy)
end


function system(definition::Logistic2Unidir)
    return DiscreteDynamicalSystem(eom_logistic2uni, definition.xi, definition)
end

Base.@kwdef struct Logistic4Chain{V, RX, RY, RZ, RW, C1, C2, C3, Σ1, Σ2, Σ3, RNG}
    xi::V = [0.1, 0.2, 0.3, 0.4]
    rx::RX = 3.9
    ry::RY = 3.6
    rz::RZ = 3.6
    rw::RW = 3.8
    c_xy::C1 = 0.4
    c_yz::C2 = 0.4
    c_zw::C3 = 0.35
    σ_xy::Σ1 = 0.05
    σ_yz::Σ2 = 0.05
    σ_zw::Σ3 = 0.05
    rng::RNG = Random.default_rng()
end

function system(definition::Logistic4Chain)
    return DiscreteDynamicalSystem(eom_logistic4_chain, definition.xi, definition)
end

function eom_logistic4_chain(u, p::Logistic4Chain, t)
    (; xi, rx, ry, rz, rw, c_xy, c_yz, c_zw, σ_xy, σ_yz, σ_zw, rng) = p
    x, y, z, w = u
    f_xy = (y +  c_xy*(x + σ_xy * rand(rng)) ) / (1 + c_xy*(1+σ_xy))
    f_yz = (z +  c_yz*(y + σ_yz * rand(rng)) ) / (1 + c_yz*(1+σ_yz))
    f_zw = (w +  c_zw*(z + σ_zw * rand(rng)) ) / (1 + c_zw*(1+σ_zw))
    dx = rx * x * (1 - x)
    dy = ry * (f_xy) * (1 - f_xy)
    dz = rz * (f_yz) * (1 - f_yz)
    dw = rw * (f_zw) * (1 - f_zw)
    return SVector{4}(dx, dy, dz, dw)
end

Base.@kwdef struct Logistic2Bidir{V, C1, C2, R1, R2, Σx, Σy, R}
    xi::V = [0.5, 0.5]
    c_xy::C1 = 0.1
    c_yx::C2 = 0.1
    r₁::R1 = 3.78
    r₂::R2 = 3.66
    σ_xy::Σx = 0.05
    σ_yx::Σy = 0.05
    rng::R = Random.default_rng()
end

function system(definition::Logistic2Bidir)
    return DiscreteDynamicalSystem(eom_logistic2bidir, definition.xi, definition)
end

# Note: Until the `eom_logistic2_bidir` function is deprecated, this function must
# be called something different; otherwise the DiscreteDynamicalSystem constructor
# doesn't work.
function eom_logistic2bidir(u, p::Logistic2Bidir, t)
    (; xi, c_xy, c_yx, r₁, r₂, σ_xy, σ_yx, rng) = p
    x, y = u
    f_xy = (y +  c_xy*(x + σ_xy * rand(rng)) ) / (1 + c_xy*(1+σ_xy))
    f_yx = (x +  c_yx*(y + σ_yx * rand(rng)) ) / (1 + c_yx*(1+σ_yx))
    dx = r₁ * (f_yx) * (1 - f_yx)
    dy = r₂ * (f_xy) * (1 - f_xy)
    return SVector{2}(dx, dy)
end