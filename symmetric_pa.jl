using Revise, CausalityTools, DynamicalSystems, TimeseriesSurrogates, Statistics, Plots

npts = 1000
c_xy = 0.5
c_yx = 0.7
sys = logistic2_bidir(c_xy = c_xy, c_yx = c_yx)
x, y = columns(trajectory(sys, npts, Ttr = 1000));
est = VisitationFrequency(RectangularBinning(3))
ηmax = 10
ηs = 1:ηmax
#x, y = rand(npts), rand(npts)

τX = estimate_delay(x, "ac_min"); dX = findmin(delay_fnn(x, τX))[2]
τY = estimate_delay(y, "ac_min"); dY = findmin(delay_fnn(y, τY))[2]

pa_xy = PredictiveAsymmetry.pa(x, y, est, ηs, τS = -τX, dS = dX, τT = -τY, dT = dY)
pa_yx = PredictiveAsymmetry.pa(y, x, est, ηs, τS = -τY, dS = dY, τT = -τX, dT = dX)

surr_type = RandomShuffle()
sx = surrogenerator(x, surr_type)
sy = surrogenerator(y, surr_type)

nsurr = 19
pa_xy_surr = zeros(ηmax, nsurr)
pa_yx_surr = zeros(ηmax, nsurr)

for i = 1:19
    pa_xy_surr[:, i] = PredictiveAsymmetry.pa(sx(), y, est, ηs, τS = -τX, dS = dX, τT = -τY, dT = dY)
    pa_xy_surr[:, i] = PredictiveAsymmetry.pa(sy(), x, est, ηs, τS = -τX, dS = dX, τT = -τY, dT = dY)
end

α = 0.05
xy_uq = [quantile(pa_xy_surr[:, i], 1 - α) for i = 1:length(ηs)]
yx_uq = [quantile(pa_yx_surr[:, i], 1 - α) for i = 1:length(ηs)]