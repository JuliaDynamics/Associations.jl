@inline @inbounds function eom_forced_lorenz_lorenz_bidir(u, p, t)
    c_xy, c_yx, c_zx, c_zy, a₁, a₂, a₃, b₁, b₂, b₃, a₃, c₁, c₂, c₃ = (p...,)
    x1, x2, x3, y1, y2, y3, z1, z2, z3 = (u...,)
    
    dx1 = a₁*(x2 - x1) + c_yx*(y1 - x1) + c_zx*(z1 - x1)
    dx2 = a₂*x1 - x1*x3 - x2
    dx3 = x1*x2 - a₃*x3
    
    dy1 = b₁*(y2 - y1) + c_yx*(x1 - y1) + c_zx*(z1 - y1)
    dy2 = b₂*y1 - y1*y3 - y2
    dy3 = y1*y2 - b₃*y3
    
    dz1 = c₁*(z2 - z1)
    dz2 = c₂*z1 - z1*z3 + z2
    dz3 = z1*z2 - c₃*z3
    
    return SVector{9}(dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3)
end

function forced_lorenz_lorenz_bidir(; u0 = rand(9), 
        c_xy = 0.1, c_yx = 0.1, # mostly synchronized for c_xy or c_yx > 0.1
        c_zx = 0.05, c_zy = 0.05, 
        a₁ = 10, a₂ = 28, a₃ = 8/3,
        b₁ = 10, b₂ = 28, b₃ = 8/3,
        c₁ = 10, c₂ = 28, c₃ = 8/3)
    ContinuousDynamicalSystem(eom_forced_lorenz_lorenz_bidir, u0, [c_xy, c_yx, c_zx, c_zy, a₁, a₂, a₃, b₁, b₂, b₃, a₃, c₁, c₂, c₃])
end

function forced_lorenz_lorenz_bidir_trajectory(npts; 
        n_transient = 2000, dt = 0.1, sample_dt = 1,
        u0 = rand(9),
        c_xy = 1.0, c_yx = 1.0, c_zx = 1.0, c_zy = 1.0, # beyond c = 2, systems syncronize
        a₁ = 10, a₂ = 28, a₃ = 8/3,
        b₁ = 10, b₂ = 28, b₃ = 8/3,
        c₁ = 10, c₂ = 28, c₃ = 8/3)

    s = forced_lorenz_lorenz_bidir(u0 = u0, 
        c_xy = c_xy, c_yx = c_yx, 
        c_zx = c_zx, c_zy = c_zy,
        a₁ = a₁, a₂ = a₂, a₃ = a₃,
        b₁ = b₁, b₂ = b₂, b₃ = b₃,
        c₁ = c₁, c₂ = c₂, c₃ = c₃)

    # the system is recorded at times t0:dt:T
    T = npts*dt*sample_dt
    o = trajectory(s, T, dt = dt, Ttr = n_transient*dt, alg = SimpleDiffEq.SimpleATsit5())[1:sample_dt:end-1, :] #alg = SimpleDiffEq.SimpleATsit5()
end

function good_forced_lorenz_lorenz_bidir_trajectory(npts; 
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
        pts = forced_lorenz_lorenz_bidir_trajectory(npts, 
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