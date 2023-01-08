export repressilator

using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: trajectory
using DynamicalSystemsBase: ContinuousDynamicalSystem

@inline @inbounds function eom_repressilator(u, p, t)
    α, α₀, n, β = p.α, p.α₀, p.n, p.β

    # pᵢ := concentration of protein repressor i
    # mᵢ := concentration of mRNA associated with pᵢ
    m₁, m₂, m₃, p₁, p₂, p₃ = u

    ṁ₁ = -m₁ + α/(1 + p₃^n) + α₀
    ṁ₂ = -m₂ + α/(1 + p₁^n) + α₀
    ṁ₃ = -m₃ + α/(1 + p₂^n) + α₀
    ṗ₁ = -β*(p₁ - m₁)
    ṗ₂ = -β*(p₂ - m₂)
    ṗ₃ = -β*(p₃ - m₃)

    return SVector{6}(ṁ₁, ṁ₂, ṁ₃, ṗ₁, ṗ₂, ṗ₃)
end

"""
    repressilator(;u₀ = rand(6), α = 10.0, α₀ = 0.0, β = 100.0,
        n = 2) → ContinuousDynamicalSystem

A six-dimensional repressilator (or repression-driven oscillator) from Elowitz & Leibler
(2000)[^Elowitz2000]. The equations are scaled to be non-dimensional.

Used in Sun & Bollt (2014) to study the performance of the causation entropy algorithm.

[^Elowitz2000]: Elowitz, M. B., & Leibler, S. (2000). A synthetic oscillatory network of
    transcriptional regulators. Nature, 403(6767), 335-338.
[^Sun2014]: Sun, J., Cafaro, C., & Bollt, E. M. (2014). Identifying the coupling structure
    in complex systems through the optimal causation entropy principle. Entropy, 16(6),
    3416-3433.
"""
function repressilator(;u₀ = rand(6), α = 10.0, α₀ = 0.0, β = 100.0, n = 2)

    p = @LArray [α, α₀, β, n] (:α, :α₀, :β, :n)
    ContinuousDynamicalSystem(eom_repressilator, u₀, p)
end
