@inline @inbounds function eom_forced_rossler_rossler_bidir(u, p, t)
    ω₁, ω₂, ω₃, c_xy, c_yx, c_zx, c_zy, a₁, a₂, a₃, b₁, b₂, b₃, a₃, c₁, c₂, c₃ = (p...,)
    x1, x2, x3, y1, y2, y3, z1, z2, z3 = (u...,)
    
    dx1 = -ω₁*(x2 + x3) + c_yx*(y1 - x1) + c_zx*(z1 - x1)
    dx2 = ω₁*x1 + a₁*x2
    dx3 = a₂ + x3*(x1 - a₃)
    
    dy1 = -ω₂*(y2 + y3) + c_xy*(x1 - y1) + c_zy*(z1 - y1)
    dy2 = ω₂*y1 + b₁*y2
    dy3 = b₂ + y3*(y1 - b₃)
    
    dz1 = -ω₂*(z2 + z3)
    dz2 = ω₂*z1 + c₁*z2
    dz3 = c₂ + z3*(z1 - c₃)
    
    return SVector{9}(dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3)
end

function forced_rossler_rossler_bidir(; u0 = rand(9), 
        ω₁ = 1.015, ω₂ = 0.985, ω₃ = 0.95,
        c_xy = 0.1, c_yx = 0.1, # mostly synchronized for c_xy or c_yx > 0.1
        c_zx = 0.05, c_zy = 0.05, 
        a₁ = 0.15, a₂ = 0.2, a₃ = 10,
        b₁ = 0.15, b₂ = 0.2, b₃ = 10,
        c₁ = 0.15, c₂ = 0.2, c₃ = 10)
    ContinuousDynamicalSystem(eom_forced_rossler_rossler_bidir, u0, [ω₁, ω₂, ω₃, c_xy, c_yx, c_zx, c_zy, a₁, a₂, a₃, b₁, b₂, b₃, a₃, c₁, c₂, c₃])
end

function forced_rossler_rossler_bidir_trajectory(npts; 
        n_transient = 2000, dt = 0.6, sample_dt = 1,
        u0 = rand(9), ω₁ = 1.015, ω₂ = 0.985, ω₃ = 0.95, 
        c_xy = 0.2, c_yx = 0.2, c_zx = 0.05, c_zy = 0.05,
        a₁ = 0.15, a₂ = 0.2, a₃ = 10,
        b₁ = 0.15, b₂ = 0.2, b₃ = 10,
        c₁ = 0.15, c₂ = 0.2, c₃ = 10)

    s = forced_rossler_rossler_bidir(u0 = u0, 
        c_xy = c_xy, c_yx = c_yx, 
        c_zx = c_zx, c_zy = c_zy,
        ω₁ = ω₁, ω₂ = ω₂, ω₃ = ω₃, 
        a₁ = a₁, a₂ = a₂, a₃ = a₃,
        b₁ = b₁, b₂ = b₂, b₃ = b₃,
        c₁ = c₁, c₂ = c₂, c₃ = c₃)

    # the system is recorded at times t0:dt:T
    T = npts*dt*sample_dt
    o = trajectory(s, T, dt = dt, Ttr = n_transient*dt, alg = SimpleDiffEq.SimpleATsit5())[1:sample_dt:end-1, :] #alg = SimpleDiffEq.SimpleATsit5()
end

function good_forced_rossler_rossler_bidir_trajectory(npts; 
        sample_dt = 1,  Ttr = 5000, dt = 0.3, # dt = 0.6 about 10 samples per period
        Da₁ = Uniform(0.12, 0.17),
        Da₂ = Uniform(0.18, 0.22),
        Da₃ = Uniform(9.0, 11.0),
        Db₁ = Uniform(0.10, 0.20),
        Db₂ = Uniform(0.18, 0.22),
        Db₃ = Uniform(9.0, 11.0),
        Dc₁ = Uniform(0.10, 0.20),
        Dc₂ = Uniform(0.18, 0.22),
        Dc₃ = Uniform(9.0, 11.0),
        Dω₁ = Uniform(0.95, 0.999), #Uniform(0.97, 1.03)
        Dω₂ = Uniform(1.001, 1.05), #Uniform(0.97, 1.03)
        Dω₃ = Uniform(0.9, 0.95), #Uniform(0.97, 1.03)

        a₁ = nothing,
        a₂ = nothing,
        a₃ = nothing,
        b₁ = nothing,
        b₂ = nothing,
        b₃ = nothing,
        c₁ = nothing,
        c₂ = nothing,
        c₃ = nothing,
        ω₁ = nothing, 
        ω₂ = nothing, 
        ω₃ = nothing, 
        c_xy = 0.2,  c_yx = 0.2,
        c_zx = 0.05, c_zy = 0.05,
        u0 = rand(9),
        n_maxtries = 300)

    n_tries = 0

    while n_tries <= n_maxtries
        ω₁ == nothing ? ω₁ = rand(Dω₁) : ω₁ = ω₁
        ω₂ == nothing ? ω₂ = rand(Dω₂) : ω₂ = ω₂
        ω₃ == nothing ? ω₃ = rand(Dω₃) : ω₃ = ω₃
        a₁ == nothing ? a₁ = rand(Da₁) : a₁ = a₁
        a₂ == nothing ? a₂ = rand(Da₂) : a₂ = a₂
        a₃ == nothing ? a₃ = rand(Da₃) : a₃ = a₃
        b₁ == nothing ? b₁ = rand(Db₁) : b₁ = b₁
        b₂ == nothing ? b₂ = rand(Db₂) : b₂ = b₂
        b₃ == nothing ? b₃ = rand(Db₃) : b₃ = b₃
        c₁ == nothing ? c₁ = rand(Dc₁) : c₁ = c₁
        c₂ == nothing ? c₂ = rand(Dc₂) : c₂ = c₂
        c₃ == nothing ? c₃ = rand(Dc₃) : c₃ = c₃
        pts = forced_rossler_rossler_bidir_trajectory(npts, 
            sample_dt = sample_dt, dt = dt, n_transient = Ttr,
            c_xy = c_xy,  c_yx = c_yx, 
            c_zx = c_zx, c_zy = c_zy,
            ω₁ = ω₁, ω₂ = ω₂, ω₃ = ω₃,
            a₁ = a₁, a₂ = a₂, a₃ = a₃,
            b₁ = b₁, b₂ = b₂, b₃ = b₃, 
            c₁ = c₁, c₂ = c₂, c₃ = c₃)
        
        if all(Matrix(pts) .< 1e10) && length(unique(pts)) > npts/2
            return pts
        end
        println("no attractor found. trying with new initial condition and parameters")
        n_tries += 1
    end
end