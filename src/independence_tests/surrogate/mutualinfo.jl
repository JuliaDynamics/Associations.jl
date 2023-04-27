function SurrogateTest(measure::MutualInformation, est::Nothing, args...; kwargs...)
    T = typeof(measure)
    txt = "Estimator not provided for measure $T. Cannot construct `SurrogateTest`\n" *
        "Do e.g. `SurrogateTest(MIShannon(), KSG2())`"
    throw(ArgumentError(txt))
end
