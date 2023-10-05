include("KLDivergence.jl")
include("joint_entropies/joint_entropies.jl")
include("HellingerDistance.jl")
include("conditional_mutual_informations/conditional_mutual_informations.jl")

abstract type MutualInformation <: BivariateInformationMeasure end
include("MIShannon.jl")
include("MITsallisMartin.jl")
include("MITsallisFuruichi.jl")
include("MIRenyiSarbu.jl")
include("MIRenyiJizba.jl")
