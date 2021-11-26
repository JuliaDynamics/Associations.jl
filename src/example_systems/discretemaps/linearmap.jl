using LabelledArrays

export linearmap1

"""
    eom_linearmap1(x, p, n) → SVector{2}

# References
Chen, Yonghong, et al. "Analyzing multiple nonlinear time series with extended
Granger causality." Physics Letters A 324.1 (2004): 26-35
"""
function eom_linearmap1(x, p, t)
    c = p[1]
    x, y = (x...,)
    t = t + 3
    dx = 3.4*x*(t - 1)*(1 - x^2*(t - 1))*exp(-x^2*(t - 1)) + 0.8*x*(t - 2) + rand(Normal(0, 0.05))
    dy = 3.4*y*(t - 1)*(1 - y^2*(t - 1))*exp(-y^2*(t - 1)) + 0.5*y*(t - 2) + c*x*(t - 2) + rand(Normal(0, 0.05))
    return SVector{2}(dx, dy)
end

"""
    linearmap1(u₀, c) → DiscreteDynamicalSystem

# References
Chen, Yonghong, et al. "Analyzing multiple nonlinear time series with extended Granger causality." Physics Letters A 324.1 (2004): 26-35
"""
function linearmap1(u₀, c)
    p = @LArray [c] (:c)
    DiscreteDynamicalSystem(eom_linearmap1, u₀, p)
end

"""
    linearmap1(;u₀ = [1, rand(2)], c = 0.5) → DiscreteDynamicalSystem

# References
Chen, Yonghong, et al. "Analyzing multiple nonlinear time series with extended Granger causality." Physics Letters A 324.1 (2004): 26-35
"""
linearmap1(;u₀ = rand(2), c = 0.5) = linearmap1(u₀, c)
