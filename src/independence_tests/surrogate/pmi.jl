function SurrogateAssociationTest(measure::PMI, est::Nothing, args...;
        kwargs...)
    T = typeof(measure)
    txt = "Estimator not provided for measure $T. Cannot construct `SurrogateAssociationTest`\n" *
        "Do e.g. `SurrogateAssociationTest(PMI(), MesnerShalizi())`"
    throw(ArgumentError(txt))
end
