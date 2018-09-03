"""
    eom_linearmap1(x, p, n) -> SVector{2}

# References
Chen, Yonghong, et al. "Analyzing multiple nonlinear time series with extended
Granger causality." Physics Letters A 324.1 (2004): 26-35
"""
function eom_linearmap1(x, p, n)
    c = p[1]
    x, y = (x...)
    dx = 3.4*x*(n - 1)*(1 - x^2*(n - 1))*exp(-x^2*(n - 1)) + 0.8*x*(n - 2)
    dy = 3.4*y*(n - 1)*(1 - y^2*(n - 1))*exp(-y^2*(n - 1)) + 0.5*y*(n - 2) + c*x*(n - 2)
    return SVector{2}(dx, dy)
end

doc"""
    linearmap1(u₀, c) -> DiscreteDynamicalSystem

"""
function linearmap1(u₀, c)
    p = [c]
    logistic_system = DiscreteDynamicalSystem(eom_linearmap1, u₀, p)
    return logistic_system
end

linearmap1(;u₀ = rand(2, c = 2.0) = linearmap1(u₀, c)
