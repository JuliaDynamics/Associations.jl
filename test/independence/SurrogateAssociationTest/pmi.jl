using Random
rng = MersenneTwister(1234)
n = 200

# Pre-discretized data
likeit = rand(rng, ["yes", "no"], n)
food = rand(rng, ["veggies", "meat", "fish"], n)
service = rand(rng, ["netflix", "hbo"], n)
est = Contingency()
nshuffles = 3

@test_throws ArgumentError SurrogateAssociationTest(PMI())

@test independence(SurrogateAssociationTest(PMI(), est; nshuffles, rng), food, likeit, service) isa SurrogateAssociationTestResult
@test independence(SurrogateAssociationTest(CMIRenyiSarbu(), est; nshuffles, rng), food, likeit, service) isa SurrogateAssociationTestResult
@test independence(SurrogateAssociationTest(CMIRenyiJizba(), est; nshuffles, rng), food, likeit, service) isa SurrogateAssociationTestResult


# Analytical tests, in the limit.
# -------------------------------
n = 10000
# Pre-discretized data
likeit = rand(rng, ["yes", "no"], n);
food = rand(rng, ["veggies", "meat", "fish"], n);
service = rand(rng, ["netflix", "hbo"], n);

α = 0.01 # pick some arbitrary significance level

# We should not be able to reject the null hypothesis `food ⫫ likeit | service`, because
# the variables are all independent.
test_cmi = independence(SurrogateAssociationTest(PMI(), est; nshuffles = 200, rng), food, likeit, service)
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
test = SurrogateAssociationTest(PMI(), est; nshuffles = 200, rng)
test_cmi = independence(test, places, experience, preferred_equipment)
@test pvalue(test_cmi) > α


nshuffles = 5
surrtest_sp = SurrogateAssociationTest(PMI(), OrdinalPatterns(); nshuffles, rng)
surrtest_vh = SurrogateAssociationTest(PMI(), ValueBinning(4); nshuffles, rng)
surrtest_dp = SurrogateAssociationTest(PMI(), Dispersion(); nshuffles, rng)

@test independence(surrtest_sp, x, y, z) isa SurrogateAssociationTestResult
@test independence(surrtest_vh, x, y, z) isa SurrogateAssociationTestResult
@test independence(surrtest_dp, x, y, z) isa SurrogateAssociationTestResult

@test independence(surrtest_sp, X, Y, Z) isa SurrogateAssociationTestResult
@test independence(surrtest_vh, X, Y, Z) isa SurrogateAssociationTestResult
@test independence(surrtest_dp, X, Y, Z) isa SurrogateAssociationTestResult
