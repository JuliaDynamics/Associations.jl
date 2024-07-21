using Random
rng = MersenneTwister(1234)
n = 100

# Pre-discretized data
likeit = rand(rng, ["yes", "no"], n)
food = rand(rng, ["veggies", "meat", "fish"], n)
service = rand(rng, ["netflix", "hbo"], n)
nshuffles = 3

@test_throws ArgumentError SurrogateAssociationTest(CMIShannon())

# Estimators
d = CodifyVariables(UniqueElements()) # discretization
est_cmi_shannon = JointProbabilities(CMIShannon(), d)
est_cmi_renyisarbu = JointProbabilities(CMIShannon(), d)
est_cmi_renyijizba = JointProbabilities(CMIShannon(), d)

# Independence tests
test_cmi_shannon = SurrogateAssociationTest(est_cmi_shannon; nshuffles, rng)
test_cmi_renyisarbu = SurrogateAssociationTest(est_cmi_renyisarbu; nshuffles, rng)
test_cmi_renyijizba = SurrogateAssociationTest(est_cmi_renyijizba; nshuffles, rng)

@test independence(test_cmi_shannon, food, likeit, service) isa SurrogateAssociationTestResult
@test independence(test_cmi_renyisarbu, food, likeit, service) isa SurrogateAssociationTestResult
@test independence(test_cmi_renyijizba, food, likeit, service) isa SurrogateAssociationTestResult

# Analytical tests, in the limit of many samples
# ----------------------------------------------
n = 1000
# Pre-discretized data
likeit = rand(rng, ["yes", "no"], n)
food = rand(rng, ["veggies", "meat", "fish"], n)
service = rand(rng, ["netflix", "hbo"], n)

α = 0.01 # pick some arbitrary significance level

# We should not be able to reject the null hypothesis `food ⫫ likeit | service`, because
# the variables are all independent.
d = CodifyVariables(UniqueElements()) # outcome space
est = JointProbabilities(CMIShannon(), d)
test = SurrogateAssociationTest(est; nshuffles = 19, rng)
test_cmi = independence(test, food, likeit, service)
@test pvalue(test_cmi) > α

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

# We should not be able to reject the null hypothesis `places ⫫ experience | preferred_equipment`, because
# places → preferred_equipment → experience, so when conditioning on the intermediate variable,
# the first and last variable in the chain should be independent.
d = CodifyVariables(UniqueElements())
est = JointProbabilities(CMIShannon(), d)
test = SurrogateAssociationTest(est; nshuffles = 19, rng)
test_cmi = independence(test, places, experience, preferred_equipment)
@test pvalue(test_cmi) > α
