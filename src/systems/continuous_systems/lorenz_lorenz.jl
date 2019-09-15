@inline @inbounds function eom_lorenz_lorenz(u, p, t)
    c_xy, c_yx, a₁, a₂, a₃, b₁, b₂, b₃ = (p...,)
    x1, x2, x3, y1, y2, y3 = (u...,)
    
    dx1 = -a₁*(x1 - x2) + c_yx*(y1 - x1)
    dx2 = -x1*x3 + a₂*x1 - x2
    dx3 = x1*x2 - a₃*x3
    dy1 = -b₁*(y1 - y2) + c_xy*(x1 - y1)
    dy2 = -y1*y3 + b₂*y1 - y2
    dy3 = y1*y2 - b₃*y3
    
    return SVector{6}(dx1, dx2, dx3, dy1, dy2, dy3)
end

function lorenz_lorenz(; u0 = rand(6), 
        c_xy = 0.2, c_yx = 0.2, 
        a₁ = 10, a₂ = 28, a₃ = 8/3, 
        b₁ = 10, b₂ = 28, b₃ = 9/3)
    ContinuousDynamicalSystem(eom_lorenz_lorenz, u0, [c_xy, c_yx, a₁, a₂, a₃, b₁, b₂, b₃])
end

function lorenzlorenz_trajectory(npts; sample_dt = 1, Ttr = 1000, dt = 0.1, 
    c_xy = 0.1, c_yx = 0.1, 
    u0 = rand(6),
    a₁ = 10, a₂ = 28, a₃ = 8/3, 
    b₁ = 10, b₂ = 28, b₃ = 9/3)

s = lorenz_lorenz(u0 = u0, c_xy = c_xy, c_yx = c_yx, a₁ = a₁, a₂ = a₂, a₃ = a₃, b₁ = b₁, b₂ = b₂, b₃ = b₃)

# the system is recorded at times t0:dt:T
T = npts*dt*sample_dt

o = trajectory(s, T, dt = dt, Ttr = Ttr*dt, alg = SimpleDiffEq.SimpleATsit5())[1:sample_dt:end-1, :]
end



# For some initial lconditions, the system wanders off and doesn't settle to an attractor. Create a function that loops until we get a good realization.
function good_lorenzlorenztrajectory(npts; 
        sample_dt = 1, 
        dt = 0.1, 
        c_xy = 0.1, 
        c_yx = 0.1,
        Da₁ = Uniform(9.5, 10.5),
        Da₂ = Uniform(27, 29),
        Da₃ = Uniform(7.5/3, 8.5/3),
        Db₁ = Uniform(9.5, 10.5),
        Db₂ = Uniform(27, 29),
        Db₃ = Uniform(7.5/3, 8.5/3),
        a₁ = nothing, 
        a₂ = nothing, 
        a₃ = nothing, 
        b₁ = nothing, 
        b₂ = nothing, 
        b₃ = nothing,
        u0 = rand(6),
        Ttr = 10000,
        n_maxtries = 300)

    n_tries = 0
    while n_tries <= n_maxtries
        a₁ == nothing ? a₁ = rand(Da₁) : nothing
        a₂ == nothing ? a₂ = rand(Da₂) : nothing
        a₃ == nothing ? a₃ = rand(Da₃) : nothing
        b₁ == nothing ? b₁ = rand(Db₁) : nothing
        b₂ == nothing ? b₂ = rand(Db₂) : nothing
        b₃ == nothing ? b₃ = rand(Db₃) : nothing
        
        pts = lorenzlorenz_trajectory(npts, 
            sample_dt = sample_dt, dt = dt, 
            c_xy = c_xy, c_yx = c_yx,
            Ttr = Ttr)
        
        M = Matrix(pts) 
        
        if all(isfinite.(M)) && all(M .< 1e10) && count(M .≈ 0) < npts*0.1 && count(abs.(M) .< 1e-10) < npts*0.1 && 
            (count(abs.(M) .< 1e-12) < npts*0.1) 
            return pts
        end
        println("no attractor found. trying with new initial condition and parameters")
        n_tries += 1
    end
end
