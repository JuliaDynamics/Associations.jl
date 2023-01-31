import TimeseriesSurrogates: surrogenerator
using TimeseriesSurrogates: RandomShuffle, SurrogateGenerator

function surrogenerator(x::AbstractDataset, rf::RandomShuffle, rng = Random.default_rng())
    n = length(x)
    idxs = collect(1:n)

    init = (
        permutation = collect(1:n),
    )

    return SurrogateGenerator(rf, x, similar(x), init, rng)
end

function (sg::SurrogateGenerator{<:RandomShuffle, <:AbstractDataset})()
    x, s, rng = sg.x, sg.s, sg.rng
    n = length(x)
    permutation = getfield.(Ref(sg.init), (:permutation))
    shuffle!(rng, permutation)
    for i in 1:n
        s[i] = x[permutation[i]]
    end
    return s
end
