export GaussianCMI
using StateSpaceSets: StateSpaceSet

"""
    GaussianCMI <: MutualInformationEstimator
    GaussianCMI(; normalize::Bool = false)

`GaussianCMI` is a parametric estimator for Shannon conditional mutual information (CMI)
[Vejmelka2008](@cite).

## Description

`GaussianCMI` estimates Shannon CMI through a sum of two mutual information terms that
each are estimated using [`GaussianMI`](@ref) (the `normalize` keyword is the same as
for [`GaussianMI`](@ref)):

```math
\\hat{I}_{Gaussian}(X; Y | Z) = \\hat{I}_{Gaussian}(X; Y, Z) - \\hat{I}_{Gaussian}(X; Z)
```
"""
Base.@kwdef struct GaussianCMI <: MutualInformationEstimator
    normalize::Bool = false
end

function estimate(measure::CMIShannon, est::GaussianCMI, x, y, z)
    YZ = StateSpaceSet(y, z)

    mi_est = GaussianMI()
    MI_x_yz = estimate(MIShannon(measure.e), mi_est, x, YZ)
    MI_x_z = estimate(MIShannon(measure.e), mi_est, x, z)

    return MI_x_yz - MI_x_z
end
