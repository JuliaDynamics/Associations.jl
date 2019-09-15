import SimpleDiffEq


""" 
    ContinuousSystemSetup(resampling = NoResampling, dt::Int = 0.01, Ttr::Int = 0, 
        sample_step::Int = 1, n_pts::Int = 100,
        diffeq = (alg = SimpleDiffEq.SimpleATsit5(), abstol = 1.0e-6, reltol = 1.0e-6),
        source, target)

Setup for causality analysis of continuous dynamical systems.

## Mandatory keyword arguments 

- **`source`**:  The variable(s) to use as `source` when calling `causality`. Usually integer(s).

- **`target`**:  The variable(s) to use as `target` when calling `causality`. Usually integer(s).

## Optional keywords

- **`resampling`**: An instance of a resampling scheme. Defaults to `NoResampling()`.

- **`dt::Int`**: The time step when iterating the system. Defaults to 0.01.

- **`Ttr::Int`**: The number of transient iterations before starting sampling. Defaults to 0.

- **`sample_step::Int`**: The sampling step in the final orbit. If `sample_step` > 1, then the 
    orbit is generated until time `T*sample_step`, after which each `sample_step`th point is 
    drawn to generate the final orbit.

- **`diffeq`**: Arguments propagated to init DifferentialEquations.jl. Defaults to
    `diffeq = (alg = SimpleDiffEq.SimpleATsit5(), abstol = 1.0e-6, reltol = 1.0e-6)`.

- **`n_pts::Int`**: The number of points in the orbit. Defaults to 100.

"""
Base.@kwdef struct ContinuousSystemSetup{R} <: DynamicalSystemSetup where R
    """ An instance of a resampling scheme (defaults to `NoResampling()`)"""
    resampling::R = NoResampling()
    
    """ The time step """
    dt = 0.01
    
    """ The number of transient iterations (in units of `dt`) before starting sampling. """
    Ttr::Int = 0
    
    """ 
    The sampling step. If `sample_step` > 1, then the orbit is generated until time 
    T*sample_step, after which each `sample_step`th point is drawn to generate the 
    final orbit.
    """
    sample_step::Int = 1  

    """ The noise level, expressed as a fraction of the standard deviation of each variable. """
    noise_level = 0

    """ The number of points in the orbit (in units of `dt`). Defaults to 100. """
    n_pts::Int = 100

    """ The variable(s) to use as `source` """
    source
    
    """ The variable(s) to use as `target` """
    target
end

export ContinuousSystemSetup