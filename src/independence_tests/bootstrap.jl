using StatsBase: sample!
using Random

function bootstrap(f::Function, x; n = 100, rng = Random.default_rng)
    s = similar(x)
    estimates = zeros(n)
    for i = 1:n
        sample!(rng, x, s; replace = true, ordered = false)
        estimates[i] = f(s)
    end
    return estimates
end
