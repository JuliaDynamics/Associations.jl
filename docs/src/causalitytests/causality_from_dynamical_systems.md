
# [Causality in dynamical systems](@id causality_dynamical_systems)

## General syntax

The syntax for quantifying directional causality between variables of dynamical systems is 

```julia
causality(system::DynamicalSystem, setup::DynamicalSystemSetup, test::CausalityTest)
```

Here,

- `system` is either a `DiscreteDynamicalSystem` or `ContinuousDynamicalSystem` instance. 
    See [`DynamicalSystems.jl`](https://github.com/JuliaDynamics/DynamicalSystems.jl) for 
    more info on how to construct and initialise dynamical systems,

- `setup` is an instance of either [`ContinuousSystemSetup`](@ref) or 
    [`DiscreteSystemSetup`](@ref). This gives the information on trajectory lengths, sampling 
    steps and solver options. It also instructs the causality estimators about which 
    variables of the dynamical system from which the time series for `source` and `target` 
    are to be sampled,

- `test` is an instance of a [causality test] (@ref causality_tests).

What happens under the hood is that the dynamical system you provided is iterated 
using the provided `setup` parameters. A trajectory for the system is then recorded and 
subsampled according to your specifications in `setup`. Then, the causality `test` is 
applied to the time series for the variables you indicated in the `setup`.

### Continuous systems

```@docs
causality(::ContinuousDynamicalSystem, ::ContinuousSystemSetup{NoResampling}, ::CausalityTest)
```

### Discrete systems

```@docs
causality(::DiscreteDynamicalSystem, ::DiscreteSystemSetup{NoResampling}, ::CausalityTest)
```

## Setting up the analysis

The setup for continuous and discrete systems is slightly different, so there are two separate 
types where you can give instructions about step length, time series length, etc, to the solvers.

### Continuous system setup

```@docs
CausalityTools.IntegrationDynamicalSystems.ContinuousSystemSetup
```

### Discrete system setup

```@docs
CausalityTools.IntegrationDynamicalSystems.DiscreteSystemSetup
```
