export GaussianCMI
using StateSpaceSets: StateSpaceSet

"""
    GaussianCMI <: MutualInformationEstimator
    GaussianCMI(definition = CMIShannon(); normalize::Bool = false)

`GaussianCMI` is a parametric [`ConditionalMutualInformationEstimator`](@ref) 
[Vejmelka2008](@cite).

## Usage

- Use with [`association`](@ref) to compute [`CMIShannon`](@ref) from input data.

## Description

`GaussianCMI` estimates Shannon CMI through a sum of two mutual information terms that
each are estimated using [`GaussianMI`](@ref) (the `normalize` keyword is the same as
for [`GaussianMI`](@ref)):

```math
\\hat{I}_{Gaussian}(X; Y | Z) = \\hat{I}_{Gaussian}(X; Y, Z) - \\hat{I}_{Gaussian}(X; Z)
```

## Example 

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 10000)
y = rand(rng, 10000) .+ x
z = rand(rng, 10000) .+ y
association(GaussianCMI(CMIShannon(base = 2)), x, z, y)
```

## Compatible definitions

- [`CMIShannon`](@ref)
"""
struct GaussianCMI{M <: ConditionalMutualInformation} <: ConditionalMutualInformationEstimator{M}
    definition::M
    normalize::Bool
end
function GaussianCMI(definition = CMIShannon(); normalize = false)
    return GaussianCMI(definition, normalize)
end

function association(est::GaussianCMI{<:CMIShannon}, x, y, z)
    YZ = StateSpaceSet(y, z)

    mi_est_modified = estimator_with_overridden_parameters(est.definition, GaussianMI())  
    MI_x_yz = association(mi_est_modified, x, YZ)
    MI_x_z = association(mi_est_modified, x, z)

    return MI_x_yz - MI_x_z
end
