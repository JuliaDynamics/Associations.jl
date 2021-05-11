using LabelledArrays

export var1

"""
    eom_var1(x, p, n) → Function

Equations of motion for a vector autoregressive system where X₁ → X₂ → X₃.
"""
function eom_var1(x, p, n)
    σ₁, σ₂, σ₃ = p[1], p[2], p[3]
    x₁, x₂, x₃ = x[1], x[2], x[3]
    θ = rand(Normal(0, σ₁))
    η = rand(Normal(0, σ₂))
    ϵ = rand(Normal(0, σ₃))

    dx₁ = θ
    dx₂ = x₁ * η
    dx₃ = 0.5*x₃ * x₂ + ϵ
    return SVector{3}(dx₁, dx₂, dx₃)
end

function var1(u₀, σ₁, σ₂, σ₃)
    p = @LArray [σ₁, σ₂, σ₃] (:σ₁, :σ₂, :σ₃)

    DiscreteDynamicalSystem(eom_var1, u₀, p)
end

"""
    var1(x, p, n) → DiscreteDynamicalSystem

Initialise a discrete vector autoregressive system where X₁ → X₂ → X₃.
"""
var1(;u₀ = rand(3), σ₁ = 1.0, σ₂ = 0.2, σ₃ = 0.3) = var1(u₀, σ₁, σ₂, σ₃)
