using Revise, CausalityTools, DynamicalSystems, TimeseriesSurrogates, Statistics, Plots, Distributions

npts = 1000
c_xy = rand(Uniform(0.2, 0.8))
c_yx = rand(Uniform(0.2, 0.8))
@show c_xy, c_yx
sys = logistic2_bidir(c_xy = c_xy, c_yx = c_yx)
x, y = columns(trajectory(sys, npts, Ttr = 1000));
ηmax = 8
ηs = 1:ηmax
#x, y = rand(npts), rand(npts)

τX = estimate_delay(x, "ac_min"); dX = findmin(delay_fnn(x, τX))[2]
τY = estimate_delay(y, "ac_min"); dY = findmin(delay_fnn(y, τY))[2]
#dTot = dX + dY + 1 # target is always 1-dimensional

nbins = floor(Int, npts^(1/(max(dX, dY) + 1)))
est = VisitationFrequency(RectangularBinning(nbins))

#pa_xy, Ȳ⁺Y⁻X⁻s_xy, Ȳ⁻Y⁺X⁻s_xy, Ȳ⁺Y⁻s_xy, Ȳ⁻Y⁺s_xy = PredictiveAsymmetry.A₁(x, y, est, ηs, τS = -τX, dS = dX, τT = -τY, dT = dY)
#pa_yx, Ȳ⁺Y⁻X⁻s_yx, Ȳ⁻Y⁺X⁻s_yx, Ȳ⁺Y⁻s_yx, Ȳ⁻Y⁺s_yx = PredictiveAsymmetry.A₁(y, x, est, ηs, τS = -τY, dS = dY, τT = -τX, dT = dX)
pa_xy = PredictiveAsymmetry.A₁(x, y, est, ηs, τS = -τX, dS = dX, τT = -τY, dT = dY)
pa_yx = PredictiveAsymmetry.A₁(y, x, est, ηs, τS = -τY, dS = dY, τT = -τX, dT = dX)


te_xy = [transferentropy(x, y, est, η𝒯 = η, τS = -τX, dS = dX, τT = -τY, dT = dY) for η in [-ηmax:-1; 1:ηmax]]
te_yx = [transferentropy(y, x, est, η𝒯 = η, τS = -τY, dS = dY, τT = -τX, dT = dX) for η in [-ηmax:-1; 1:ηmax]]

surr_type = RandomShuffle()#BlockShuffle(floor(Int, npts / 8))
sx = surrogenerator(x, surr_type)
sy = surrogenerator(y, surr_type)

nsurr = 100
pa_xy_surr = zeros(ηmax, nsurr)
pa_yx_surr = zeros(ηmax, nsurr)

# pe = plot(xlabel = "η", ylabel = "Entropy (bits)", legend = :right)
# plot!(ηs, Ȳ⁺Y⁻X⁻s_xy .-  Ȳ⁺Y⁻s_xy, label = "X → Y | H(X⁻|Ȳ⁺Y⁻)", ls = :solid, c = :black)
# plot!(ηs, Ȳ⁻Y⁺X⁻s_xy .-  Ȳ⁻Y⁺s_xy, label = "X → Y | H(X⁻|Ȳ⁻Y⁺)", ls = :dash, c = :black)

# plot!(ηs, Ȳ⁺Y⁻X⁻s_yx .-  Ȳ⁺Y⁻s_yx, label = "Y → X | H(X⁻|Ȳ⁺Y⁻)", ls = :solid, c = :red)
# plot!(ηs, Ȳ⁻Y⁺X⁻s_yx .-  Ȳ⁻Y⁺s_yx, label = "Y → Y | H(X⁻|Ȳ⁻Y⁺)", ls = :dash, c = :red)

# plot!(ηs, Ȳ⁺Y⁻X⁻s_xy, label = "X → Y | Ȳ⁺Y⁻X⁻", ls = :dash, c = :black)
# plot!(ηs, Ȳ⁻Y⁺X⁻s_xy, label = "X → Y | Ȳ⁻Y⁺X⁻", ls = :dot, c = :black)
# plot!(ηs, Ȳ⁺Y⁻s_xy, label = "X → Y | Ȳ⁺Y⁻", ls = :dashdot, c = :black)
# plot!(ηs, Ȳ⁻Y⁺s_xy, label = "X → Y | Ȳ⁻Y⁺", ls = :solid, c = :black)

# plot!(ηs, Ȳ⁺Y⁻X⁻s_yx, label = "Y → X | Ȳ⁺Y⁻X⁻", ls = :dash, c = :red)
# plot!(ηs, Ȳ⁻Y⁺X⁻s_yx, label = "Y → X | Ȳ⁻Y⁺X⁻", ls = :dot, c = :red)
# plot!(ηs, Ȳ⁺Y⁻s_yx, label = "Y → X | Ȳ⁺Y⁻", ls = :dashdot, c = :red)
# plot!(ηs, Ȳ⁻Y⁺s_yx, label = "Y → X | Ȳ⁻Y⁺", ls = :solid, c = :red)

for i = 1:nsurr
    pa_xy_surr[:, i] = PredictiveAsymmetry.A₁(sx(), y, est, ηs, τS = -τX, dS = dX, τT = -τY, dT = dY)
    pa_yx_surr[:, i] = PredictiveAsymmetry.A₁(sy(), x, est, ηs, τS = -τY, dS = dY, τT = -τX, dT = dX)
end

α = 0.01
uq_xy = [quantile(pa_xy_surr[:, i], 1 - α) for i = 1:length(ηs)]
uq_yx = [quantile(pa_yx_surr[:, i], 1 - α) for i = 1:length(ηs)]

ymax = maximum(abs.([pa_xy; pa_yx; uq_xy; uq_yx]))*1.1
pa = plot(xlabel = "η", ylabel = "PA", ylims = (-ymax, ymax))
plot!(ηs, pa_xy, label = "X → Y", c = :black)
plot!(ηs, uq_xy, label = "", c = :black, ls = :dash)

plot!(ηs, pa_yx, label = "Y → X", c = :blue)
plot!(ηs, uq_yx, label = "", c = :blue, ls = :dash)
hline!([0], label = "", ls = :dot, c = :grey)

px = plot(xlabel = "Time step", ylabel = "Value")
plot!(x, label = "x", c = :black)
py = plot(xlabel = "Time step", ylabel = "Value")
plot!(y, label = "y", c = :blue)
pts = plot(px, py, layout = grid(2, 1))

pte = plot(xlabel = "η", ylabel = "TE")
plot!([-ηmax:-1; 1:ηmax], te_xy, c = :black, label = "X → Y")
plot!([-ηmax:-1; 1:ηmax], te_yx, c = :blue, label = "Y → X")
vline!([0], label = "", ls = :dash, c = :grey)
plot(pts, plot(pte, pa, layout = grid(2, 1)), layout = grid(1, 2), size = (1200, 800))