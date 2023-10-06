abstract type ConditionalEntropy <: MultivariateInformationMeasure end

min_inputs_vars(::ConditionalEntropy) = 2
max_inputs_vars(::ConditionalEntropy) = 2

include("CEShannon.jl")
include("CETsallisAbe.jl")
include("CETsallisFuruichi.jl")