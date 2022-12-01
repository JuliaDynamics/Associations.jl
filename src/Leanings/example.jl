using CausalityTools, DynamicalSystems, Plots, Statistics, Distributions; gr()

npts, Ttr = 10000, 5000
x, y = columns(trajectory(henon2(c_xy = 0.0), npts - 1, Ttr = Ttr))
xr, yr = rand(Uniform(-1, 1), npts), rand(Uniform(-1, 1), npts)
plot()
plot!(x, label = "x")
plot!(y, label = "y")

import Entropies: symbolize, SymbolicPermutation

m = 2; τ = -1
n = 5
weighted = false
nreps = 50
ls = 0:10
leans_wn = zeros(nreps, length(ls))
leans_henon_c0 = zeros(nreps, length(ls))
leans_henon_c2 = zeros(nreps, length(ls))

function partition(x, n::Int)
    xmin = minimum(x) # the reference point
    Δx = (maximum(x) - xmin) / n
    [floor(Int, (xᵢ - xmin) / Δx) for xᵢ in x]
end

for i = 1:nreps 
    xr, yr = rand(Uniform(-1, 1), npts), rand(Uniform(-1, 1), npts)
    X = partition(xr, n)#symbolize(xr, SymbolicPermutation(m = m, τ = τ))
    Y = partition(yr, n)#symbolize(yr, SymbolicPermutation(m = m, τ = τ))
    for (j, l) in enumerate(ls)
        leans_wn[i, j] = lean(X, Y, l, weighted = weighted)
    end

    x, y = columns(trajectory(henon2(c_xy = 0.0), npts - 1, Ttr = Ttr))
    X = partition(x, n)#symbolize(x, SymbolicPermutation(m = m, τ = τ))
    Y = partition(y, n)#symbolize(y, SymbolicPermutation(m = m, τ = τ))
    for (j, l) in enumerate(ls)
        leans_henon_c0[i, j] = lean(X, Y, l, weighted = weighted)
    end

    x, y = columns(trajectory(henon2(c_xy = 2.0), npts - 1, Ttr = Ttr))
    X = partition(x, n)#symbolize(x, SymbolicPermutation(m = m, τ = τ))
    Y = partition(y, n)#symbolize(y, SymbolicPermutation(m = m, τ = τ))
    for (j, l) in enumerate(ls)
        leans_henon_c2[i, j] = lean(X, Y, l, weighted = weighted)
    end
end

mean_leans_wn = dropdims(mean(leans_wn, dims = 1), dims = 1)
mean_leans_henon_c0 = dropdims(mean(leans_henon_c0, dims = 1), dims = 1)
mean_leans_henon_c2 = dropdims(mean(leans_henon_c2, dims = 1), dims = 1)

plot(xlabel = "l", ylabel = "λx→y", ylims = (-1.1, 1.1))
#plot!(ls, mean_leans_wn, label = "white noise", marker = :star)
plot!(ls, mean_leans_henon_c0, label = "henon (c_xy = 0.0)", marker = :circle)
#plot!(ls, mean_leans_henon_c2, label = "henon (c_xy = 2.0)", marker = :star5)