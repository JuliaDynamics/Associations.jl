
using Test
using Random
rng = MersenneTwister(1234)
n = 100

# Pre-discretized data
likeit = rand(rng, ["yes", "no"], n)
food = rand(rng, ["veggies", "meat", "fish"], n)
service = rand(rng, ["netflix", "hbo"], n)
nshuffles = 3

@test_throws ArgumentError SurrogateAssociationTest(PartialMutualInformation())

# estimator 
d = CodifyVariables(UniqueElements())
est_pmi = JointProbabilities(PartialMutualInformation(), d)

# Tests 
test = SurrogateAssociationTest(est_pmi; nshuffles, rng)
@test independence(test, food, likeit, service) isa SurrogateAssociationTestResult


# Analytical tests, in the limit.
# -------------------------------
n = 3000
# Pre-discretized data
likeit = rand(rng, ["yes", "no"], n);
food = rand(rng, ["veggies", "meat", "fish"], n);
service = rand(rng, ["netflix", "hbo"], n);

α = 0.01 # pick some arbitrary significance level

# We should not be able to reject the null hypothesis `food ⫫ likeit | service`, because
# the variables are all independent.
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
test_pmi = independence(test, places, experience, preferred_equipment)
@test pvalue(test_cmi) > α


# Numeric tests
x, y, z = rand(rng, n), rand(rng, n), rand(rng, n)
X, Y, Z = StateSpaceSet(x), StateSpaceSet(y), StateSpaceSet(z)
nshuffles = 19
d_ord = CodifyVariables(OrdinalPatterns())
d_disp = CodifyVariables(Dispersion())
d_bin = CodifyVariables(ValueBinning(4))
est_ord = JointProbabilities(PartialMutualInformation(), d_ord)
est_disp = JointProbabilities(PartialMutualInformation(), d_disp)
est_bin = JointProbabilities(PartialMutualInformation(), d_bin)

surrtest_ord = SurrogateAssociationTest(est_ord; nshuffles, rng)
surrtest_disp = SurrogateAssociationTest(est_disp; nshuffles, rng)
surrtest_bin = SurrogateAssociationTest(est_bin; nshuffles, rng)

@test independence(surrtest_ord, x, y, z) isa SurrogateAssociationTestResult
@test independence(surrtest_disp, x, y, z) isa SurrogateAssociationTestResult
@test independence(surrtest_bin, x, y, z) isa SurrogateAssociationTestResult

@test independence(surrtest_ord, X, Y, Z) isa SurrogateAssociationTestResult
@test independence(surrtest_disp, X, Y, Z) isa SurrogateAssociationTestResult
@test independence(surrtest_bin, X, Y, Z) isa SurrogateAssociationTestResult
