using Random
rng = MersenneTwister(1234)
n = 100

# Pre-discretized data
likeit = rand(rng, ["yes", "no"], n)
food = rand(rng, ["veggies", "meat", "fish"], n)
service = rand(rng, ["netflix", "hbo"], n)
est = 
nshuffles = 3
d = CodifyVariables(UniqueElements())
test_mi_s = SurrogateAssociationTest(JointProbabilities(MIShannon(), d); nshuffles, rng)
test_mi_rj = SurrogateAssociationTest(JointProbabilities(MIRenyiJizba(), d); nshuffles, rng)
test_mi_rs = SurrogateAssociationTest(JointProbabilities(MIRenyiSarbu(), d); nshuffles, rng)
test_mi_tf = SurrogateAssociationTest(JointProbabilities(MITsallisFuruichi(), d); nshuffles, rng)
test_mi_tm = SurrogateAssociationTest(JointProbabilities(MITsallisMartin(), d); nshuffles, rng)

@test independence(test_mi_s, food, likeit) isa SurrogateAssociationTestResult
@test independence(test_mi_rj, food, likeit) isa SurrogateAssociationTestResult
@test independence(test_mi_rs, food, likeit) isa SurrogateAssociationTestResult
@test independence(test_mi_tf, food, likeit) isa SurrogateAssociationTestResult
@test independence(test_mi_tm, food, likeit) isa SurrogateAssociationTestResult

# Analytical tests, in the limit.
# -------------------------------
n = 100000
α = 0.02 # pick some arbitrary significance level

# Simulate a survey where the place a person grew up controls how many times they
# fell while going skiing. The control happens through an intermediate variable
# `preferred_equipment`, which indicates what type of physical activity the
# person has engaged with. For this example, we should be able to reject
# places ⫫ experience, but not reject places ⫫ experience | preferred_equipment
places = rand(rng, ["city", "countryside", "under a rock"], n);
preferred_equipment = map(places) do place
    if cmp(place, "city") == 1
        return rand(rng, ["skateboard", "bmx bike"])
    elseif cmp(place, "countryside") == 1
        return rand(rng, ["sled", "snowcarpet"])
    else
        return rand(rng, ["private jet", "car"])
    end
end;
experience = map(preferred_equipment) do equipment
    if equipment ∈ ["skateboard", "bmx bike"]
        return "didn't fall"
    elseif equipment ∈ ["sled", "snowcarpet"]
        return "fell 3 times or less"
    else
        return "fell uncontably many times"
    end
end;

# We should be able to reject the null hypothesis of `places ⫫ experience`.
test_mi = independence(SurrogateAssociationTest(MIShannon(), est; nshuffles, rng), places, experience)
@test pvalue(test_mi) < α
