include("cov.jl")
include("kde.jl")
include("cca.jl")
include("multidimensional_surrogates.jl")
include("extensions.jl")
include("transformations.jl")

function logq0(q)
    if q == 1.0
        return x -> zero(x)
    else
        return x -> (x^(1 - q) - 1)/(1 - q)
    end
end
