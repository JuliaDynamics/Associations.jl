
""" 
    DiscreteSystemSetup(resampling = NoResampling, dt::Int = 1, Ttr::Int = 0, 
        n_pts = 100, source, target)

Setup for causality analysis of discrete dynamical systems.

## Mandatory keyword arguments 

- **`source`**:  The variable(s) to use as `source` when calling `causality`. Usually integer(s).

- **`target`**:  The variable(s) to use as `target` when calling `causality`. Usually integer(s).

## Optional keywords

- **`resampling`**: An instance of a resampling scheme. Defaults to `NoResampling()`.

- **`dt::Int`**: The time step when iterating the system.

- **`Ttr::Int`**: The number of transient iterations before starting sampling.

- **`diffeq`**: Arguments propagated to init DifferentialEquations.jl. Defaults to
    `diffeq = (alg = SimpleDiffEq.SimpleATsit5(), abstol = 1.0e-6, reltol = 1.0e-6)`.

- **`n_pts::Int`**: The number of points in the orbit. 
"""
Base.@kwdef struct DiscreteSystemSetup{R} <: DynamicalSystemSetup where R
    """ An instance of a resampling scheme (defaults to `NoResampling()`)"""
    resampling::R = NoResampling()
    
    """ The time step """
    dt::Int = 1
    
    """ The number of transient iterations before starting sampling """
    Ttr::Int = 0
    
    """ The number of points in the orbit """
    n_pts = 100
    
    """ The variable(s) to use as `source` """
    source
    
    """ The variable(s) to use as `target` """
    target
end

export DiscreteSystemSetup