using DynamicalSystems 

# Observational noise 
x = ones(10, 2) .+ rand(Normal(0, 0.1), 10, 2)
@test add_observational_noise!(x, 10) isa typeof(x)
@test add_observational_noise!(Dataset(x), 10) isa Dataset

#####################
# Continuous systems
#####################

# Rössler-Lorenz
# --------------
@test RosslerLorenzUnidir() isa RosslerLorenzUnidir
@test ContinuousDynamicalSystem(RosslerLorenzUnidir()) isa ContinuousDynamicalSystem
@test rand(RosslerLorenzUnidir) isa RosslerLorenzUnidir
@test trajectory(RosslerLorenzUnidir(), 100) isa Dataset