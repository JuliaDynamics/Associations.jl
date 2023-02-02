# [Independence testing](@id quickstart_independence)

## Mutual information (categorical)

```@example
using CausalityTools
# Simulate 
n = 1000
preference = rand(["yes", "no"], n)
food = rand(["veggies", "meat", "fish"], n)
test = SurrogateTest(MIShannon(), Contingency())
independence(test, preference, food)
```

## Conditional mutual information (categorical)

Here, we simulate a survey at a ski resort. The data are such that the place a person
grew up is associated with how many times they fell while going skiing. The control
happens through an intermediate variable `preferred_equipment`, which indicates what
type of physical activity the person has engaged with in the past. Some activities
like skateboarding leads to better overall balance, so people that are good on
a skateboard also don't fall, and people that to less challenging activities fall
more often.

We should be able to reject `places ⫫ experience`, but not reject
`places ⫫ experience | preferred_equipment`.  Let's see if we can detect these
relationships using (conditional) mutual information.

```@docs
using CausalityTools
n = 100000
α = 0.02 # pick some arbitrary significance level

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

test_mi = independence(SurrogateTest(MIShannon(), est; nsurr), places, experience)
test_cmi = independence(SurrogateTest(CMIShannon(), est; nsurr), places, experience, preferred_equipment)
```
