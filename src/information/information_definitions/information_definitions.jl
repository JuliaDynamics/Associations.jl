include("divergences_and_distances/divergences_and_distances.jl")
include("joint_entropies/joint_entropies.jl")
include("mutual_informations/mutual_informations.jl")
include("conditional_mutual_informations/conditional_mutual_informations.jl")

abstract type MutualInformation <: BivariateInformationMeasure end
include("MIShannon.jl")
include("MITsallisMartin.jl")
include("MITsallisFuruichi.jl")
include("MIRenyiSarbu.jl")
include("MIRenyiJizba.jl")
