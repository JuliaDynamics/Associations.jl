include("cov.jl")
include("kde.jl")
include("cca.jl")

log0(base, x) = x == 0 ? 0 : log(base, x)
