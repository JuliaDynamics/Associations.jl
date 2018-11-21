"""
    eom_chuacircuit_nscroll_sine(u, p, t) -> SVector{3}

Equations of motion for n-scroll chaotic attractors from an adjusted Chua
system [1].

# References
1. Tang, Wallace KS, et al. "Generation of n-scroll attractors via
sine function." IEEE Transactions on Circuits and Systems I:
Fundamental Theory and Applications 48.11 (2001): 1369-1372.

"""
function eom_chuacircuit_nscroll_sine(u, p, t)
    α, β, γ, a, b, c, σx, σy, σz = (p...,)
    x, y, z = (u...,)

    n::Int = c + 1
    if x >= 2*a*c
        fx = (b*pi/2*a)*(x - 2*a*c)
    elseif -2*a*c < x < 2*a*c
        d = ifelse(isodd(n), pi, 0)
        fx = -b*sin((pi*x/2*a) + d)
    elseif x <= -2*a*c
        fx = (b*pi/2*a)*(x + 2*a*c)
    end

    ηx = σx == 0 ? 0 : rand(Normal(0, σx))
    ηy = σy == 0 ? 0 : rand(Normal(0, σy))
    ηz = σz == 0 ? 0 : rand(Normal(0, σz))

    dx = α*(y - fx) + ηx
    dy = x - y + z + ηy
    dz = -β*y - γ*z + ηz
    return SVector{3}(dx, dy, dz)
end

"""
    chua_nscroll_sinefunc(u₀, α, β, γ, a, b, c::Int)

Generate n-scroll chaotic attractors from an adjusted Chua system [1].

# References
1. Tang, Wallace KS, et al. "Generation of n-scroll attractors via
sine function." IEEE Transactions on Circuits and Systems I:
Fundamental Theory and Applications 48.11 (2001): 1369-1372.
"""
function chuacircuit_nscroll_sine(u₀, α, β, γ, a, b, c::Int, σx, σy, σz)
    p = [α, β, γ, a, b, c, σx, σy, σz]
    ContinuousDynamicalSystem(eom_chuacircuit_nscroll_sine, u₀, p)
end
chuacircuit_nscroll_sine(;u₀ = [0.0, 0.0, 0.28695],
        α = 10.814, β = 14, γ = 0, a = 1.3, b = 0.11, c = 2,
        σx = 0.0, σy = 0.0, σz = 0.0) =
    chuacircuit_nscroll_sine(u₀, α, β, γ, a, b, c, σx, σy, σz)
