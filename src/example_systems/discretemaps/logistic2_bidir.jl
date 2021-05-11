using LabelledArrays

export logistic2_bidir

"""
    logistic2_bidir(u₀, c_xy, c_yx, r₁, r₂, σ_xy, σ_yx)

Equations of motion for a bidirectional logistic model for the chaotic 
population dynamics of two interacting species. This system is from [1], 
and is given by 

```math
\\begin{align}
x(t+1) &= r_1 f_{yx}^{t}(1 - f_{yx}^{t}) \\
y(t+1) &= r_2 f_{xy}^{t}(1 - f_{xy}^{t}) \\
f_{xy}^t &= \\dfrac{y(t) + c_{xy}(x(t) + \\sigma_{xy} \\xi_{xy}^t )}{1 + c_{xy} (1 + \\sigma_{xy} )} \\ 
f_{yx}^t &= \\dfrac{x(t) + c_{yx}(y(t) + \\sigma_{yx} \\xi_{yx}^t )}{1 + c_{yx} (1 + \\sigma_{yx} )},
\\end{align}
```

where the coupling strength ``c_{xy}`` controls how strongly species ``x`` influences species 
``y``, and vice versa for ``c_{yx}``. To simulate time-varying influence of unobserved 
processes, we use the dynamical noise terms ``\\xi_{xy}^t`` and ``\\xi_{yx}^t``, drawn from a 
uniform distribution with support on ``[0, 1]``. If ``\\sigma_{xy} > 0``, then the influence 
of ``x`` on ``y`` is masked by dynamical noise equivalent to ``\\sigma_{xy} \\xi_{xy}^{t}`` at 
the ``t``-th iteration of the map, and vice versa for ``\\sigma_{yx}``.

## References 

1. Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation 
    using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
"""
function eom_logistic2_bidir(dx, x, p, n)
    
    # c_xy is the coupling from x to y
    # c_yx is the coupling from y to x
    # σ_yx is the dynamical noise from y to x
    # σ_xy is the dynamical noise from y to x
    c_xy, c_yx, r₁, r₂, σ_xy, σ_yx = (p...,)
    
    ξ₁ = rand() # random number from flat distribution on [0, 1]
    ξ₂ = rand() # random number from flat distribution on [0, 1]
    x, y = x[1], x[2]
    
    f_xy = (y +  c_xy*(x + σ_xy*ξ₁) ) / (1 + c_xy*(1+σ_xy))
    f_yx = (x +  c_yx*(y + σ_yx*ξ₂) ) / (1 + c_yx*(1+σ_yx))

    dx[1] = r₁ * (f_yx) * (1 - f_yx)
    dx[2] = r₂ * (f_xy) * (1 - f_xy)
    return
end

function logistic2_bidir(u₀, c_xy, c_yx, r₁, r₂, σ_xy, σ_yx)
    p = @LArray [c_xy, c_yx, r₁, r₂, σ_xy, σ_yx] (:c_xy, :c_yx, :r₁, :r₂, :σ_xy, :σ_yx)
    DiscreteDynamicalSystem(eom_logistic2_bidir, u₀, p)
end

"""
    logistic2_bidir(;u₀ = rand(2), c_xy = 0.1, c_yx = 0.1, 
        r₁ = 3.78, r₂ = 3.66, σ_xy = 0.05, σ_yx = 0.05)

A bidirectional logistic model for the chaotic population dynamics of two interacting 
species [1].

## Equations of motion 

The equations of motion are 

```math
\\begin{align}
x(t+1) &= r_1 f_{yx}^{t}(1 - f_{yx}^{t}) \\\\
y(t+1) &= r_2 f_{xy}^{t}(1 - f_{xy}^{t}) \\\\
f_{xy}^t &= \\dfrac{y(t) + c_{xy}(x(t) + \\sigma_{xy} \\xi_{xy}^t )}{1 + c_{xy} (1 + \\sigma_{xy} )} \\\\ 
f_{yx}^t &= \\dfrac{x(t) + c_{yx}(y(t) + \\sigma_{yx} \\xi_{yx}^t )}{1 + c_{yx} (1 + \\sigma_{yx} )},
\\end{align}
```

where the coupling strength ``c_{xy}`` controls how strongly species ``x`` influences species 
``y``, and vice versa for ``c_{yx}``. To simulate time-varying influence of unobserved 
processes, we use the dynamical noise terms ``\\xi_{xy}^t`` and ``\\xi_{yx}^t``, drawn from a 
uniform distribution with support on ``[0, 1]``. If ``\\sigma_{xy} > 0``, then the influence 
of ``x`` on ``y`` is masked by dynamical noise equivalent to ``\\sigma_{xy} \\xi_{xy}^{t}`` at 
the ``t``-th iteration of the map, and vice versa for ``\\sigma_{yx}``.

## References 

1. Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation 
    using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
"""
logistic2_bidir(;u₀ = rand(2), c_xy = 0.1, c_yx = 0.1, 
    r₁ = 3.78, r₂ = 3.66, σ_xy = 0.05, σ_yx = 0.05) =
    logistic2_bidir(u₀, c_xy, c_yx, r₁, r₂, σ_xy, σ_yx)
