@inline @inbounds function eom_rossler_rossler_bidir(u, p, t)
    ω₁, ω₂, c_xy, c_yx, a₁, a₂, a₃, b₁, b₂, b₃ = (p...,)
    x1, x2, x3, y1, y2, y3 = (u...,)
    
    dx1 = -ω₁*(x2 + x3) + c_yx*(y1 - x1)
    dx2 = ω₁*x1 + a₁*x2
    dx3 = a₂ + x3*(x1 - a₃)
    
    dy1 = -ω₂*(y2 + y3) + c_xy*(x1 - y1)
    dy2 = ω₂*y1 + b₁*y2
    dy3 = b₂ + y3*(y1 - b₃)
    
    return SVector{6}(dx1, dx2, dx3, dy1, dy2, dy3)
end

function rossler_rossler_bidir(; u0 = rand(6), 
        ω₁ = 1.015, ω₂ = 0.985, 
        c_xy = 0.1, c_yx = 0.1, # mostly synchronized for c_xy or c_yx > 0.1
        a₁ = 0.15, a₂ = 0.2, a₃ = 10,
        b₁ = 0.15, b₂ = 0.2, b₃ = 10)
    ContinuousDynamicalSystem(eom_rossler_rossler_bidir, u0, [ω₁, ω₂, c_xy, c_yx, a₁, a₂, a₃, b₁, b₂, b₃])
end


function rossler_rossler_bidir_trajectory(npts, sample_dt; n_transient = 2000, dt = 0.2,
        u0 = rand(6), ω₁ = 1.015, ω₂ = 0.985, 
        c_xy = 0.2, c_yx = 0.2,
        a₁ = 0.15, a₂ = 0.2, a₃ = 10,
        b₁ = 0.15, b₂ = 0.2, b₃ = 10)

    s = rossler_rossler_bidir(u0 = u0, 
        c_xy = c_xy, c_yx = c_yx,
        ω₁ = ω₁, ω₂ = ω₂, 
        a₁ = a₁, a₂ = a₂, a₃ = a₃,
        b₁ = b₁, b₂ = b₂, b₃ = b₃)

    # the system is recorded at times t0:dt:T
    T = npts*dt*sample_dt
    o = trajectory(s, T, dt = dt, Ttr = n_transient*dt, alg = SimpleDiffEq.SimpleATsit5())[1:sample_dt:end-1, :] #alg = SimpleDiffEq.SimpleATsit5()
end

function good_rossler_rossler_bidir_trajectory(npts; sample_dt = 1,
        Da₁ = Uniform(0.12, 0.17),
        Da₂ = Uniform(0.18, 0.22),
        Da₃ = Uniform(9.0, 11.0),
        Db₁ = Uniform(0.10, 0.20),
        Db₂ = Uniform(0.18, 0.22),
        Db₃ = Uniform(9.0, 11.0),
        Dω₁ = Uniform(0.95, 0.999), #Uniform(0.97, 1.03)
        Dω₂ = Uniform(1.001, 1.05), #Uniform(0.97, 1.03)
        a₁ = nothing,
        a₂ = nothing,
        a₃ = nothing,
        b₁ = nothing,
        b₂ = nothing,
        b₃ = nothing,
        ω₁ = nothing, 
        ω₂ = nothing, 
        c_xy = 0.2,  c_yx = 0.2,
        Ttr = 5000, dt = 0.2, 
        u0 = rand(6),
        n_maxtries = 300)

    n_tries = 0

    while n_tries <= n_maxtries
        ω₁ == nothing ? ω₁ = rand(Dω₁) : ω₁ = ω₁
        ω₂ == nothing ? ω₂ = rand(Dω₂) : ω₂ = ω₂
        a₁ == nothing ? a₁ = rand(Da₁) : a₁ = a₁
        a₂ == nothing ? a₂ = rand(Da₂) : a₂ = a₂
        a₃ == nothing ? a₃ = rand(Da₃) : a₃ = a₃
        b₁ == nothing ? b₁ = rand(Db₁) : b₁ = b₁
        b₂ == nothing ? b₂ = rand(Db₂) : b₂ = b₂
        b₃ == nothing ? b₃ = rand(Db₃) : b₃ = b₃
        
        pts = rossler_rossler_bidir_trajectory(npts, sample_dt, dt = dt, n_transient = Ttr,
            c_xy = c_xy,  c_yx = c_yx, ω₁ = ω₁, ω₂ = ω₂, a₁ = a₁, a₂ = a₂,a₃ = a₃,b₁ = b₁,  b₂ = b₂, b₃ = b₃)
        
        if all(Matrix(pts) .< 1e10) && length(unique(pts)) > npts/2
            return pts
        end
        println("no attractor found. trying with new initial condition and parameters")
        n_tries += 1
    end
end