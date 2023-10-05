abstract type MutualInformation <: BivariateInformationMeasure end

include("MIShannon.jl")
include("MITsallisMartin.jl")
include("MITsallisFuruichi.jl")
include("MIRenyiSarbu.jl")
include("MIRenyiJizba.jl")