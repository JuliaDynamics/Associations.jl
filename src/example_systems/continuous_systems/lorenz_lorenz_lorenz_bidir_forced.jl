using LabelledArrays

export lorenz_lorenz_lorenz_bidir_forced

@inline @inbounds function eom_lorenz_lorenz_lorenz_bidir_forced(u, p, t)
    c_xy, c_yx, c_zx, c_zy, a₁, a₂, a₃, b₁, b₂, b₃, a₃, c₁, c₂, c₃ = (p...,)
    x₁, x₂, x₃, y₁, y₂, y₃, z₁, z₂, z₃ = (u...,)
   
    dx₁ = -a₁*(x₁ - x₂) + c_yx*(y₁ - x₁) + c_zx*(z₁ - x₁)
    dx₂ = -x₁*x₃ + a₂*x₁ - x₂
    dx₃ = x₁*x₂ - a₃*x₃
    
    dy₁ = -b₁*(y₁ - y₂) + c_xy*(x₁ - y₁) + c_zy*(z₁ - y₁)
    dy₂ = -y₁*y₃ + b₂*y₁ - y₂
    dy₃ = y₁*y₂ - b₃*y₃
    
    dz₁ = -c₁*(z₁ - z₂)
    dz₂ = -z₁*z₃ + c₂*z₁ - z₂
    dz₃ = z₁*z₂ - c₃*z₃
    
    return SVector{9}(dx₁, dx₂, dx₃, dy₁, dy₁, dy₃, dz₁, dz₂, dz₃)
end

""" 
    lorenz_lorenz_lorenz_bidir_forced(; u0 = rand(9), 
        c_xy = 0.1, c_yx = 0.1,
        c_zx = 0.05, c_zy = 0.05, 
        a₁ = 10, a₂ = 28, a₃ = 8/3,
        b₁ = 10, b₂ = 28, b₃ = 8/3,
        c₁ = 10, c₂ = 28, c₃ = 8/3)

Initialise a system consisting of two bidirectionally coupled 3D Lorenz 
systems forced by an external 3D Lorenz system, giving a 9D system.

## Equations of motion 

The dynamics is generated by the following vector field

```math 
\\begin{aligned}
\\dot{x_1} &= - a_1 (x_1 - x_2) + c_{yx}(y_1 - x_1) + c_{zx}(z_1 - x_1) \\\\
\\dot{x_2} &= - x_1 x_3 + a_2 x_1 - x_2 \\\\
\\dot{x_3} &= x_1 x_2 - a_3 x_3 \\\\
\\dot{y_1} &= -b_1 (y_1 - y_2) + c_{xy} (x_1 - y_1) + c_{zy}(z_1 - y_1) \\\\
\\dot{y_2} &= - y_1 y_3 + b_2 y_1 - y_2 \\\\
\\dot{y_3} &= y_1 y_2 - b_3 y_3 \\\\
\\dot{z_1} &= - c_1 (z_1 - z_2) \\\\
\\dot{z_2} &= - z_1 z_3 + c_2 z_1 - z_2 \\\\
\\dot{z_3} &= z_1 z_2 - c_3 z_3 
\\end{aligned}
```

## References 

1. Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from 
    multivariate flows by the joint distance distribution." Chaos: An 
    Interdisciplinary Journal of Nonlinear Science 28.7 (2018): 075302.
"""
function lorenz_lorenz_lorenz_bidir_forced(; u0 = rand(9), 
        c_xy = 0.1, c_yx = 0.1,
        c_zx = 0.05, c_zy = 0.05, 
        a₁ = 10, a₂ = 28, a₃ = 8/3,
        b₁ = 10, b₂ = 28, b₃ = 8/3,
        c₁ = 10, c₂ = 28, c₃ = 8/3)

    p = @LArray [c_xy, c_yx, c_zx, c_zy, a₁, a₂, a₃, b₁, b₂, b₃, a₃, c₁, c₂, c₃] (:c_xy, :c_yx, :c_zx, :c_zy, :a₁, :a₂, :a₃, :b₁, :b₂, :b₃, :a₃, :c₁, :c₂, :c₃)
    ContinuousDynamicalSystem(eom_lorenz_lorenz_lorenz_bidir_forced, u0, p)
end

function lorenz_lorenz_lorenz_bidir_forced_trajectory(npts; 
        n_transient = 2000, dt = 0.1, sample_dt = 1,
        u0 = rand(9),
        c_xy = 1.0, c_yx = 1.0, c_zx = 1.0, c_zy = 1.0, # beyond c = 2, systems syncronize
        a₁ = 10, a₂ = 28, a₃ = 8/3,
        b₁ = 10, b₂ = 28, b₃ = 8/3,
        c₁ = 10, c₂ = 28, c₃ = 8/3)

    s = lorenz_lorenz_lorenz_bidir_forced(u0 = u0, 
        c_xy = c_xy, c_yx = c_yx, 
        c_zx = c_zx, c_zy = c_zy,
        a₁ = a₁, a₂ = a₂, a₃ = a₃,
        b₁ = b₁, b₂ = b₂, b₃ = b₃,
        c₁ = c₁, c₂ = c₂, c₃ = c₃)

    # the system is recorded at times t0:dt:T
    T = npts*dt*sample_dt
    o = trajectory(s, T, dt = dt, Ttr = n_transient*dt, alg = SimpleDiffEq.SimpleATsit5())[1:sample_dt:end-1, :] #alg = SimpleDiffEq.SimpleATsit5()
end

function good_lorenz_lorenz_lorenz_bidir_forced_trajectory(npts; 
        sample_dt = 1,  Ttr = 5000, dt = 0.1, 
        Da₁ = Uniform(9.5, 10.5),
        Da₂ = Uniform(27.5, 28.5),
        Da₃ = Uniform(7.5/3, 8.5/3),
        Db₁ = Uniform(9.5, 10.5),
        Db₂ = Uniform(27.5, 28.5),
        Db₃ = Uniform(7.5/3, 8.5/3),
        Dc₁ = Uniform(9.5, 10.5),
        Dc₂ = Uniform(27.5, 28.5),
        Dc₃ = Uniform(7.5/3, 8.5/3),

        a₁ = nothing,
        a₂ = nothing,
        a₃ = nothing,
        b₁ = nothing,
        b₂ = nothing,
        b₃ = nothing,
        c₁ = nothing,
        c₂ = nothing,
        c₃ = nothing,
        c_xy = 0.2,  c_yx = 0.2,
        c_zx = 0.05, c_zy = 0.05,
        u0 = [rand(Uniform(0, 10)) for i = 1:9],
        n_maxtries = 300)

    n_tries = 0

    while n_tries <= n_maxtries
        a₁ == nothing ? a₁ = rand(Da₁) : a₁ = a₁
        a₂ == nothing ? a₂ = rand(Da₂) : a₂ = a₂
        a₃ == nothing ? a₃ = rand(Da₃) : a₃ = a₃
        b₁ == nothing ? b₁ = rand(Db₁) : b₁ = b₁
        b₂ == nothing ? b₂ = rand(Db₂) : b₂ = b₂
        b₃ == nothing ? b₃ = rand(Db₃) : b₃ = b₃
        c₁ == nothing ? c₁ = rand(Dc₁) : c₁ = c₁
        c₂ == nothing ? c₂ = rand(Dc₂) : c₂ = c₂
        c₃ == nothing ? c₃ = rand(Dc₃) : c₃ = c₃
        pts = lorenz_lorenz_lorenz_bidir_forced(npts, 
            sample_dt = sample_dt, dt = dt, n_transient = Ttr,
            c_xy = c_xy,  c_yx = c_yx, 
            c_zx = c_zx, c_zy = c_zy,
            a₁ = a₁, a₂ = a₂, a₃ = a₃,
            b₁ = b₁, b₂ = b₂, b₃ = b₃, 
            c₁ = c₁, c₂ = c₂, c₃ = c₃)
        
        
        
        if all(Matrix(pts) .< 1e9) #&& length(unique(pts)) < length(pts)*0.8
            return pts
        end
        println("no attractor found. trying with new initial condition and parameters")
        n_tries += 1
    end
end