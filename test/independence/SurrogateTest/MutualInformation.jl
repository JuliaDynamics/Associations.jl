
n = 100

# Pre-discretized data
likeit = rand(["yes", "no"], n)
food = rand(["veggies", "meat", "fish"], n)
service = rand(["netflix", "hbo"], n)
est = Contingency()
nshuffles = 3

@test_throws ArgumentError SurrogateTest(MIShannon())

@test independence(SurrogateTest(MIShannon(), est; nshuffles), food, likeit) isa SurrogateTestResult
@test independence(SurrogateTest(MIRenyiJizba(), est; nshuffles), food, likeit) isa SurrogateTestResult
@test independence(SurrogateTest(MIRenyiSarbu(), est; nshuffles), food, likeit) isa SurrogateTestResult
@test independence(SurrogateTest(MITsallisFuruichi(), est; nshuffles), food, likeit) isa SurrogateTestResult
@test independence(SurrogateTest(MITsallisMartin(), est; nshuffles), food, likeit) isa SurrogateTestResult

# Analytical tests, in the limit.
# -------------------------------
n = 100000
α = 0.02 # pick some arbitrary significance level

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

# We should be able to reject the null hypothesis of `places ⫫ experience`.
test_mi = independence(SurrogateTest(MIShannon(), est; nshuffles), places, experience)
@test pvalue(test_mi) < α
