using LabelledArrays

export logistic2_unidir

"""
    eom_logistic2(dx, x, p, n) → function

Equations of motions for a system consisting of two coupled logistic maps where
X unidirectionally influences Y [1].


## Equations of motion

The equations of motion are 

```math
\\begin{aligned}
x(t+1) &= r_1 x(t)(1 - x(t)) \\\\
y(t+1) &= r_2 f(x,y)(1 - f(x,y)),
\\end{aligned}
```

with

```math
\\begin{aligned}
f(x,y) = \\dfrac{y + \\frac{c_{xy}(x \\xi )}{2}}{1 + \\frac{c_{xy}}{2}(1+ \\sigma )}
\\end{aligned}
```

The parameter `c_xy` controls how strong the dynamical forcing is. If `σ > 0`,
dynamical noise masking the influence of  `x` on `y` equivalent to
``\\sigma \\cdot \\xi`` is added at each iteration. Here,``\\xi`` is a draw from a
flat distribution on ``[0, 1]``. Thus, setting `σ = 0.05` is equivalent to
add dynamical noise corresponding to a maximum of ``5 \\%`` of the possible
range of values of the logistic map.

1. Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation 
    using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
"""

function eom_logistic2_unidir(dx, x, p, n)
    c_xy, r₁, r₂, σ = (p...,)
    ξ = rand() # random number from flat distribution on [0, 1]
    x, y = x[1], x[2]
    f_xy = (y +  (c_xy*(x + σ*ξ)/2) ) / (1 + (c_xy/2)*(1+σ))

    dx[1] = r₁ * x * (1 - x)
    dx[2] = r₂ * (f_xy) * (1 - f_xy)
    return
end

function logistic2_unidir(u₀, c_xy, r₁, r₂, σ)
    p = @LArray [c_xy, r₁, r₂, σ] (:c_xy, :r₁, :r₂, :σ)
    DiscreteDynamicalSystem(eom_logistic2_unidir, u₀, p)
end

"""
    logistic2(;u₀ = rand(2), c_xy = 0.1, σ = 0.05,
        r₁ = 3.78, r₂ = 3.66) → DiscreteDynamicalSystem

Initialise a system consisting of two coupled logistic maps where X
unidirectionally influences Y. By default, the parameters `r₁` and `r₂` are set
to values yielding chaotic behaviour.

## Equations of motion

The equations of motion are

```math
\\begin{aligned}
x(t+1) &= r_1 x(t)(1 - x(t)) \\\\
y(t+1) &= r_2 f(x,y)(1 - f(x,y)),
\\end{aligned}
```

with

```math
\\begin{aligned}
f(x,y) = \\dfrac{y + \\frac{c_{xy}(x \\xi )}{2}}{1 + \\frac{c_{xy}}{2}(1+ \\sigma )}
\\end{aligned}
```

The parameter `c_xy` controls how strong the dynamical forcing is. If `σ > 0`,
dynamical noise masking the influence of  `x` on `y` equivalent to
``\\sigma \\cdot \\xi`` is added at each iteration. Here,``\\xi`` is a draw from a
flat distribution on ``[0, 1]``. Thus, setting `σ = 0.05` is equivalent to
add dynamical noise corresponding to a maximum of ``5 \\%`` of the possible
range of values of the logistic map.

## References

1. Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation 
    using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
"""
logistic2_unidir(;u₀ = rand(2), c_xy = 0.1, r₁ = 3.78, r₂ = 3.66, σ = 0.05) =
    logistic2_unidir(u₀, c_xy, r₁, r₂, σ)

#To get chaotic realisation, check that the orbit doesn't settle to a few unique values
function good_logistic_unidir_trajectory(npts::Int; 
        Ttr = 1000, dt = 1,
        c_xy = 0.5,
        Dr₁ = Uniform(3.6, 4.0), 
        Dr₂ = Uniform(3.6, 4.0), 
        σ = 0.0, 
        n_maxtries = 300)
    
    n_tries = 0
    while n_tries <= n_maxtries
        s = logistic2_unidir(u₀ = rand(2),
            c_xy = c_xy, 
            σ = σ,
            r₁ = rand(Dr₁), 
            r₂ = rand(Dr₂))
        
        o = trajectory(s, npts * dt - 1, Ttr = Ttr, dt = dt)
        
        # Ensure there are not too many repeated values, so we don't have trivial behaviour
        
        if length(unique(o[:, 1])) > npts * 0.9 && length(unique(o[:, 2])) > npts * 0.9 
            return o
        end
        
        n_tries += 1
    end
end
