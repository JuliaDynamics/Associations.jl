abstract type ConditionalMutualInformation <: MultivariateInformationMeasure end

min_inputs_vars(::ConditionalMutualInformation) = 3
max_inputs_vars(::ConditionalMutualInformation) = 3

include("CMIShannon.jl")
include("CMITsallis.jl")
include("CMIRenyiJizba.jl")
include("CMIRenyiPoczos.jl")
include("CMIRenyiSarbu.jl")