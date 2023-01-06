include("cov.jl")
include("kde.jl")
include("cca.jl")

import ComplexityMeasures: log_with_base

function log_with_base(base, x)
    if x == zero(x)
        return x -> zero(x)
    else
        log_with_base(base)
    end
end

function logq0(q)
    if q == 1.0
        return x -> zero(x)
    else
        return x -> (x^(1 - q) - 1)/(1 - q)
    end
end
