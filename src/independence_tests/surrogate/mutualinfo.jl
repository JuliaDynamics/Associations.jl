function SurrogateAssociationTest(measure::MutualInformation, est::Nothing, args...; kwargs...)
    T = typeof(measure)
    txt = "Estimator not provided for measure $T. Cannot construct `SurrogateAssociationTest`\n" *
        "Do e.g. `SurrogateAssociationTest(MIShannon(), KSG2())`"
    throw(ArgumentError(txt))
end
