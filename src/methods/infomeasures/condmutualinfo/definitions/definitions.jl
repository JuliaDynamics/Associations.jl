# There are multiple ways of estimating CMI.
# `ConditionalMutualInformationEstimator`s estimates a specific CMI directly, using some
# specialized method, which is documented in their docstrings.
#
# It is also possible to compute CMI by decomposing it as a sum of less
# complicated quantities such as entropy, mutual information or some divergence.
# For each type of CMI (e.g. Shannon CMI), there are multiple estimation
# routes: a sum of four entropies, a sum of two mutual informations, and more.
#
# Additionally, there are multiple *definitions* of the same type of CMI, that in the
# literature are referred to by the same name. For example, Tsallis CMI can be defined in
# more than one way. Because definitions differ, we can't simply make
# *one* function that estimates all mutual information. I've therefore made the
# design choice that `condmutualinfo` dispatches on `CMIDefinition`s to indicate which
# formula/definition is being computed.
#
# This works really well with meta-functions like `independence`, because we can just
# plug in any `CMIDefinition` subtype with any `CMI` subtype, provide an appropriate
# estimator, and hence measure conditional independence between variables `X`, `Y`, `Z`.

# ----------------------------------------------------------------
# Shannon CMI (`CMIShannon`)
# ----------------------------------------------------------------
# There are multiple ways of expressing the Shannon CMI and estimating it.
# Each struct below defines a unique way of doing so, and implements relevant methods
# for `estimate`, if not using generic dispatch as defined in `common_dispatch.jl`.
include("CMIDefinitionShannonMI2.jl")
include("CMIDefinitionShannonH4.jl")

# ----------------------------------------------------------------
# Tsallis CMI (`CMITsallis`)
# ----------------------------------------------------------------
# There are multiple ways of expressing the Tsallis CMI and estimating it.
# Each struct below defines a unique way of doing so, and implements relevant methods
# for `estimate`, if not using generic dispatch as defined in `common_dispatch.jl`.

# ----------------------------------------------------------------
# Renyi CMI
# ----------------------------------------------------------------
# There are multiple ways of expressing the Rénryi mutual information and estimating it.
# Each struct below defines a unique way of doing so, and implements relevant methods
# for `estimate`, if not using generic dispatch as defined in `common_dispatch.jl`.
include("CMIDefinitionRenyiSarbu.jl")
