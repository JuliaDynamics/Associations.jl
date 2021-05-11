using LabelledArrays

export logistic4

"""
    eom_logistic4(u, p, t)

Equations of motion for a 4D transitive causal chain of unidirectionally
coupled logistic maps, where ``y_1 \\to y_2 \\to y_3 \\to y_4`` [1].

*Note: With the default parameters which are as in [1], for some initial conditions, 
this system wanders off to ``\\pm \\infty`` for some of the variables. Make sure that 
you have a good realisation before using the orbit for anything.**

## Equations of motion

```math 
\\begin{aligned}
y_1(t+1) &= y_1(t)(r_1 - r_1 y_1) \\\\
y_2(t+1) &= y_2(t)(r_2 - c_{12} y_1 - r_2 y_2) \\\\
y_3(t+1) &= y_3(t)(r_3 - c_{23} y_2 - r_3 y_3) \\\\
y_4(t+1) &= y_4(t)(r_4 - c_{34} y_3 - r_4 y_4)
\\end{aligned}
```

## References

1. Ye, Hao, et al. "Distinguishing time-delayed causal interactions using
    convergent cross mapping." Scientific reports 5 (2015): 14750
"""
function eom_logistic4(u, p, t)
    r₁, r₂, r₃, r₄, c₁₂, c₂₃, c₃₄  = (p...,)
    y₁, y₂, y₃, y₄ = (u...,)

    dy₁ = y₁*(r₁ - r₁*y₁)
    dy₂ = y₂*(r₂ - c₁₂*y₁ - r₂*y₂)
    dy₃ = y₃*(r₃ - c₂₃*y₂ - r₃*y₃)
    dy₄ = y₄*(r₄ - c₃₄*y₃ - r₄*y₄)
    return SVector{4}(dy₁, dy₂, dy₃, dy₄)
end


function logistic4(u₀, r₁, r₂, r₃, r₄, c₁₂, c₂₃, c₃₄)
    p = @LArray [r₁, r₂, r₃, r₄, c₁₂, c₂₃, c₃₄] (:r₁, :r₂, :r₃, :r₄, :c₁₂, :c₂₃, :c₃₄)
    DiscreteDynamicalSystem(eom_logistic4, u₀, p)
end

"""
    logistic4(;u₀ = rand(4), r₁ = 3.9, r₂ = 3.6, r₃ = 3.6, r₄ = 3.8,
        c₁₂ = 0.4, c₂₃ = 0.4, c₃₄ = 0.35) → DiscreteDynamicalSystem

Initialise a system of a transitive chain of four unidirectionally coupled
logistic maps, where ``y_1 \\to y_2 \\to y_3 \\to y_4`` [1]. Default 
parameters are as in [1].

*Note: With the default parameters which are as in [1], for some initial conditions, 
this system wanders off to ``\\pm \\infty`` for some of the variables. Make sure that 
you have a good realisation before using the orbit for anything.*

## Equations of motion

```math 
\\begin{aligned}
y_1(t+1) &= y_1(t)(r_1 - r_1 y_1) \\\\
y_2(t+1) &= y_2(t)(r_2 - c_{12} y_1 - r_2 y_2) \\\\
y_3(t+1) &= y_3(t)(r_3 - c_{23} y_2 - r_3 y_3) \\\\
y_4(t+1) &= y_4(t)(r_4 - c_{34} y_3 - r_4 y_4)
\\end{aligned}
```

## References

1. Ye, Hao, et al. "Distinguishing time-delayed causal interactions using
    convergent cross mapping." Scientific reports 5 (2015): 14750
"""
logistic4(;u₀ = rand(4),
            r₁ = 3.9, r₂ = 3.6, r₃ = 3.6, r₄ = 3.8,
            c₁₂ = 0.4, c₂₃ = 0.4, c₃₄ = 0.35) =
    logistic4(u₀, r₁, r₂, r₃, r₄, c₁₂, c₂₃, c₃₄)
