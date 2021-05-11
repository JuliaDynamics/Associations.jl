using LabelledArrays

export chuacircuits_driven

"""
    eom_chuacircuits_driven(u, p, t) → SVector{6}

Equations of motion for a bidirectionally coupled system consisting of two
driven Chua circuits [1].

## References

1. Murali, K., and M. Lakshmanan. "Chaotic dynamics of the driven Chua's 
    circuit." IEEE Transactions on Circuits and Systems I Fundamental
    Theory and Applications 40.11 (1993): 836-840.
"""
function eom_chuacircuits_driven(u, p, t)
    α₁, α₂, β₁, β₂, F₁, F₂, ω₁, ω₂, ϵ₁, ϵ₂, m₀, m₁, σ = (p...,)
    x₁, y₁, z₁ = (u[1:3]...,)
    x₂, y₂, z₂ = (u[4:6]...,)

    # Dynamical noise
    if σ == 0
        ξ = 0
    else
        ξ = rand(Normal(0, σ))
    end

    hx₁ = m₁*x₁ + 0.5*(m₀ - m₁)*(abs(x₁+1) - abs(x₁-1))
    hx₂ = m₁*x₂ + 0.5*(m₀ - m₁)*(abs(x₂+1) - abs(x₂-1))

    dx₁ = α₁*(y₁-hx₁) - α₁*ϵ₂*(x₁ - x₂)
    dy₁ = x₁-y₁+z₁
    dz₁ = -β₁*y₁ + F₁*sin(ω₁*t) + ξ

    dx₂ = α₂*(y₂-hx₂) - α₂*ϵ₁*(x₁ - x₂)
    dy₂ = x₂-y₂+z₂
    dz₂ = -β₂*y₂ + F₂*sin(ω₂*t) + ξ
    SVector{6}(dx₁, dy₁, dz₁, dx₂, dy₂, dz₂)
end

"""
    chuacircuits_driven(u₀, α₁, α₂, β₁, β₂, F₁, F₂,
                ω₁, ω₂, ϵ₁, ϵ₂, m₀, m₁, σ) → ContinuousDynamicalSystem

Initialize a bidirectionally coupled system consisting of two driven Chua
circuits [1], X₁ and X₂. The subsystems are mutually coupled by a linear
resistor, where `ϵ₁` controls the influence of X₁ on X₂, and `ϵ₂` controls the
influence of X₂ on X₁. The parameters for the subsystems are
set equal to each other, as in the original paper, but can be tuned
individually for each subsystem.

## References

1. Murali, K., and M. Lakshmanan. "Chaotic dynamics of the driven Chua's 
    circuit." IEEE Transactions on Circuits and Systems I Fundamental 
    Theory and Applications 40.11 (1993): 836-840.
"""
function chuacircuits_driven(u₀, α₁, α₂, β₁, β₂, F₁, F₂,
                                      ω₁, ω₂, ϵ₁, ϵ₂, m₀, m₁, σ)
    p = @LArray [α₁, α₂, β₁, β₂, F₁, F₂, ω₁, ω₂, ϵ₁, ϵ₂, m₀, m₁, σ] (:α₁, :α₂, :β₁, :β₂, :F₁, :F₂, :ω₁, :ω₂, :ϵ₁, :ϵ₂, :m₀, :m₁, :σ)
    ContinuousDynamicalSystem(eom_chuacircuits_driven, u₀, p)
end

"""
    chuacircuits_driven(;u₀ = [0.1, 0.1, 0.2, 0.15, 0.15, 0.22],
        α₁ = 7.0, α₂ = 7.0, β₁ = 14.286, β₂ = 14.286,
        F₁ = 1.5, F₂ = 1.5, ω₁ = 3.0, ω₂ = 3.0,
        σ = 0.1, ϵ₁ = 0.1, ϵ₂ = 0.1, m₀ = -1/7, m₁ = 2/7) → ContinuousDynamicalSystem

Initialize a bidirectionally coupled system consisting of two driven Chua
circuits [1], X₁ and X₂. The subsystems are mutually coupled by a linear
resistor, where `ϵ₁` controls the influence of X₁ on X₂, and `ϵ₂` controls the
influence of X₂ on X₁. The parameters for the subsystems are
set equal to each other, as in the original paper, but can be tuned
individually for each subsystem.

## References

1. Murali, K., and M. Lakshmanan. "Chaotic dynamics of the driven Chua's 
    circuit." IEEE Transactions on Circuits and Systems I Fundamental 
    Theory and Applications 40.11 (1993): 836-840.
"""
chuacircuits_driven(;u₀ = [0.1, 0.1, 0.2, 0.15, 0.15, 0.22],
                            α₁ = 7.0, α₂ = 7.0,
                            β₁ = 14.286, β₂ = 14.286,
                            F₁ = 1.5, F₂ = 1.5,
                            ω₁ = 3.0, ω₂ = 3.0,
                            σ = 0.1,
                            ϵ₁ = 0.1, ϵ₂ = 0.1,
                            m₀ = -1/7, m₁ = 2/7) =
    chuacircuits_driven(u₀, α₁, α₂, β₁, β₂, F₁, F₂, ω₁, ω₂, ϵ₁, ϵ₂, m₀, m₁, σ)
