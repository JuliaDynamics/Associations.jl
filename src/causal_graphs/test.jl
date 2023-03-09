
# using DelayEmbeddings
# npts = 800
# Δt = 0.1

npts = 1000
#X = columns(trajectory(system(Logistic4Chain(xi = rand(4))), npts, Ttr = 10000));
Δt = 0.1
X = columns(trajectory(system(
  RosslerBidir6(xi = rand(6))), npts * Δt; Δt, Ttr = 10000));
#X = columns(trajectory(system(RosslerLorenzUnidir6(xi = rand(6), c_xy = 1.5)), npts * Δt, Δt = Δt, Ttr = 10000));
w = [estimate_delay(x, "mi_min", 1:30) for x in X] |> maximum
@show w
alg = OPA(τmax = 0, α = 0.05, m = 5, est = FPVP(k = 1, w = w))
a = infer_graph(alg, X, verbose = true)
