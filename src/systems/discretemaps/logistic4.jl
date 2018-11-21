
"""
    eom_logistic4(u, p, t)

Equations of motion for a 4D transitive causal chain of unidirectionally
coupled logistic maps, where X1 -> X2 -> X3 -> X4 [1].

# References
Ye, Hao, et al. "Distinguishing time-delayed causal interactions using
convergent cross mapping." Scientific reports 5 (2015): 14750
"""
function eom_logistic4(u, p, t)
    r₁, r₂, r₃, r₄, c₂, c₃, c₄  = (p...,)
    y₁, y₂, y₃, y₄ = (u...,)

    dy₁ = y₁*(r₁ - r₁*y₁)
    dy₂ = y₂*(r₂ - c₂*y₁ - r₂*y₂)
    dy₃ = y₃*(r₃ - c₃*y₂ - r₃*y₃)
    dy₄ = y₄*(r₄ - c₄*y₃ - r₄*y₄)
    return SVector{4}(dy₁, dy₂, dy₃, dy₄)
end


function logistic4(u₀, r₁, r₂, r₃, r₄, c₂, c₃, c₄)
    p = [r₁, r₂, r₃, r₄, c₂, c₃, c₄]
    DiscreteDynamicalSystem(eom_logistic4, u₀, p)
end

"""
    logistic4(;u₀ = rand(4), r₁ = 3.9, r₂ = 3.6, r₃ = 3.6, r₄ = 3.8,
        c₂ = 0.4, c₃ = 0.4, c₄ = 0.35) -> DiscreteDynamicalSystem

Initialise a system of a a transitive chain of four unidirectionally coupled
logistic maps, where X1 -> X2 -> X3 -> X4 [1]. Default parameters are as in
Ye et al. (2015).

# References
Ye, Hao, et al. "Distinguishing time-delayed causal interactions using convergent cross mapping." Scientific reports 5 (2015): 14750
"""
logistic4(;u₀ = rand(4),
            r₁ = 3.9, r₂ = 3.6, r₃ = 3.6, r₄ = 3.8,
            c₂ = 0.4, c₃ = 0.4, c₄ = 0.35) =
    logistic4(u₀, r₁, r₂, r₃, r₄, c₂, c₃, c₄)
