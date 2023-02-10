
n = 200

# Pre-discretized data
likeit = rand(["yes", "no"], n)
food = rand(["veggies", "meat", "fish"], n)
service = rand(["netflix", "hbo"], n)
est = Contingency()
nshuffles = 3

@test_throws ArgumentError SurrogateTest(CMIShannon())

@test independence(SurrogateTest(CMIShannon(), est; nshuffles), food, likeit, service) isa SurrogateTestResult
@test independence(SurrogateTest(CMIRenyiSarbu(), est; nshuffles), food, likeit, service) isa SurrogateTestResult
@test independence(SurrogateTest(CMIRenyiJizba(), est; nshuffles), food, likeit, service) isa SurrogateTestResult


# Analytical tests, in the limit.
# -------------------------------
n = 10000

# Pre-discretized data
likeit = rand(["yes", "no"], n)
food = rand(["veggies", "meat", "fish"], n)
service = rand(["netflix", "hbo"], n)

α = 0.01 # pick some arbitrary significance level

# We should not be able to reject the null hypothesis `food ⫫ likeit | service`, because
# the variables are all independent.
test_cmi = independence(SurrogateTest(CMIShannon(), est; nshuffles = 200), food, likeit, service)
@test pvalue(test_cmi) > α

# Simulate a survey where the place a person grew up controls how many times they
# fell while going skiing. The control happens through an intermediate variable
# `preferred_equipment`, which indicates what type of physical activity the
# person has engaged with. For this example, we should be able to reject
# places ⫫ experience, but not reject places ⫫ experience | preferred_equipment

places = rand(["city", "countryside", "under a rock"], n);
preferred_equipment = map(places) do place
    if cmp(place, "city") == 1
        return rand(["skateboard", "bmx bike"])
    elseif cmp(place, "countryside") == 1
        return rand(["sled", "snowcarpet"])
    else
        return rand(["private jet", "car"])
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
test = SurrogateTest(CMIShannon(), est; nshuffles = 200)
test_cmi = independence(test, places, experience, preferred_equipment)
@test pvalue(test_cmi) > α
