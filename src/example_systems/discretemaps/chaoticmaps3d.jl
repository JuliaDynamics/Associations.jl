using LabelledArrays

export chaoticmaps

# """
#     eom_linearmap1(x, p, n) → SVector{2}

# # References
# Chen, Yonghong, et al. "Analyzing multiple nonlinear time series with extended
# Granger causality." Physics Letters A 324.1 (2004): 26-35
# """
# function eom_linearmap1(x, p, t)
#     c = p[1]
#     x, y = (x...,)
#     dx = 3.4*x*(t - 1)* (1 - x^2*(t - 1)) * exp((-x^2)*(t - 1)) + 0.8*x*(t - 2)
#     dy = 3.4*y*(t - 1)* (1 - y^2*(t - 1)) * exp((-y^2)*(t - 1)) + 0.5*y*(t - 2) + c*x*(t - 2)
#     return SVector{2}(dx, dy)
# end

# """
#     linearmap1(u₀, c) → DiscreteDynamicalSystem

# # References
# Chen, Yonghong, et al. "Analyzing multiple nonlinear time series with extended Granger causality." Physics Letters A 324.1 (2004): 26-35
# """
# function linearmap1(u₀, c)
#     p = @LArray [c] (:c)
#     logistic_system = DiscreteDynamicalSystem(eom_linearmap1, u₀, p)
#     return logistic_system
# end

# """
#     linearmap1(;u₀ = [1, rand(2)], c = 2.0) → DiscreteDynamicalSystem

# # References
# Chen, Yonghong, et al. "Analyzing multiple nonlinear time series with extended Granger causality." Physics Letters A 324.1 (2004): 26-35
# """
# linearmap1(;u₀ = rand(2), c = 0.5) = linearmap1(u₀, c)



"""
    eom_chaoticmaps(x, p, n) → SVector{3}

# References
Chen, Yonghong, et al. "Analyzing multiple nonlinear time series with extended
Granger causality." Physics Letters A 324.1 (2004): 26-35
"""
function eom_chaoticmaps(x, p, t)
    c = p[1]
    x, y, z = (x...)
    t = t + 2
    dx = 3.4*x*(t - 1) * (1 - x^2*(t - 1)) * exp((-x^2)*(t - 1))
    dy = 3.4*y*(t - 1) * (1 - y^2*(t - 1)) * exp((-y^2)*(t - 1)) + 0.5*x*(t - 1)
    dz = 3.4*z*(t - 1) * (1 - z^2*(t - 1)) * exp((-z^2)*(t - 1)) + c*x*(t - 1) + 0.3*y*(t - 1)

    return SVector{3}(dx, dy, dz)
end

function chaoticmaps(u₀, c)
    p = [c]
    cmap = DiscreteDynamicalSystem(eom_chaoticmaps, u₀, p)
    return cmap
end

"""
    chaoticmaps(;u₀ = rand(3), c = 2.0) → DiscreteDynamicalSystem

# References
Chen, Yonghong, et al. "Analyzing multiple nonlinear time series with extended Granger causality." Physics Letters A 324.1 (2004): 26-35
"""
chaoticmaps(;u₀ = rand(3), c = 0.5) = chaoticmaps(u₀, c)
