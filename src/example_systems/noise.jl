
"""
    noise_uu(n::Int, lo = - 1, hi = 1)

Generate a signal consisting of `n` steps of uncorrelated uniform noise from 
a uniform distribution on `[lo, hi]`.
"""
function noise_uu(n::Int; lo = - 1, hi = 1)
    u = Uniform(-lo, hi)
    rand(u, n)
end


"""
    noise_ug(n::Int; μ = 0, σ = 1)

Generate a signal consisting of `n` steps of uncorrelated Gaussian noise from
a normal distribution with mean `μ` and standard deviation `σ`.
"""
function noise_ug(n::Int; μ = 0, σ = 1)
    d = Normal(μ, σ)
    rand(d, n)
end

"""
    noise_brownian(n::Int; lo = - 1, hi = 1)
    noise_brownian(d::Distribution, n::Int)

Generate a signal consisting of `n` steps of Brownian noise, generated as
the zero-mean and unit standard deviation normalised cumulative sum of noise
generated from a uniform distribution on `[lo, hi]`. Optionally, a distribution
`d` from which to sample can be provided.

## Examples

```julia
# Based on uncorrelated uniform noise
noise_brownian(100)
noise_brownian(100, lo = -2, hi = 2)
noise_brownian(Uniform(-3, 3), 100)

# Based on uncorrelated Gaussian noise
μ, σ = 0, 2
noise_brownian(Normal(μ, σ), 100)
```
"""
function noise_brownian(n::Int; lo = - 1, hi = 1)
    u = Uniform(lo, hi)
    xs = cumsum(rand(u, n))
    (xs .- mean(xs)) ./ std(xs)
end

function noise_brownian(d::Distribution, n::Int)
    xs = cumsum(rand(d, n))
    (xs .- mean(xs)) ./ (std(xs))
end

export noise_uu, noise_ug, noise_brownian