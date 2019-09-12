import StatsBase

x, y = rand(300), rand(300)

####################
# Single binnings
####################
single_binning = RectangularBinning(0.5)

@test causality(x, y, VisitationFrequencyTest(binning = single_binning, ηs = 5)) isa Array{Float64,0}
@test causality(x, y, VisitationFrequencyTest(binning = single_binning, ηs = 2:10)) isa Vector

####################################################################################################
# Multiple binnings (summarised to a single value per η using the `binning_summary_statistic`).
####################################################################################################
multiple_binnings = [RectangularBinning(n) for n in 3:5]

causality(x, y, VisitationFrequencyTest(binning = multiple_binnings, 
    binning_summary_statistic = StatsBase.mean, ηs = 1)) isa Array{Float64,0}

causality(x, y, VisitationFrequencyTest(binning = multiple_binnings, 
    binning_summary_statistic = StatsBase.mean, ηs = 2:10)) isa Vector


################################################################
# Integration with UncertainData.jl
################################################################
uvals_x = [UncertainValue(Normal, rand(Normal(0, 5)), abs(rand(Normal(0, 3)))) for i = 1:100]
uvals_y = [UncertainValue(Normal, rand(Normal(0, 5)), abs(rand(Normal(0, 3)))) for i = 1:100];

# UncertainValueDataset
UVX = UncertainValueDataset(uvals_x)
UVY = UncertainValueDataset(uvals_y)

# UncertainIndexDataset
UVX_idx = UncertainIndexDataset(uvals_x)
UVY_idx = UncertainIndexDataset(uvals_y)

# Real-valued vectors
x = resample.(uvals_x);
y = resample.(uvals_y);
binnings = [RectangularBinning(i) for i = 3:6]
vftest = VisitationFrequencyTest(binning = binnings, ηs = 1)

@test causality(x, y, vftest) isa Array{<:Real, 0}
@test causality(uvals_x, uvals_y, vftest) isa Array{<:Real, 0}
@test causality(x, uvals_y, vftest) isa Array{<:Real, 0}
@test causality(uvals_x, y, vftest) isa Array{<:Real, 0}
@test causality(UVX, UVY, vftest) isa Array{<:Real, 0}
@test causality(x, UVY, vftest) isa Array{<:Real, 0}
@test causality(UVX, y, vftest) isa Array{<:Real, 0}

vftest = VisitationFrequencyTest(binning = binnings, ηs = 1:5)

@test causality(x, y, vftest) isa Array{T, 1} where T
@test causality(uvals_x, uvals_y, vftest) isa Array{T, 1} where T
@test causality(x, uvals_y, vftest) isa Array{T, 1} where T
@test causality(uvals_x, y, vftest) isa Array{T, 1} where T
@test causality(UVX, UVY, vftest) isa Array{T, 1} where T
@test causality(x, UVY, vftest) isa Array{T, 1} where T
@test causality(UVX, y, vftest) isa Array{T, 1} where T

@test causality(x, y, vftest) |> length == 5
@test causality(uvals_x, uvals_y, vftest) |> length == 5
@test causality(x, uvals_y, vftest) |> length == 5
@test causality(uvals_x, y, vftest) |> length == 5
@test causality(UVX, UVY, vftest) |> length == 5
@test causality(x, UVY, vftest) |> length == 5
@test causality(UVX, y, vftest) |> length == 5