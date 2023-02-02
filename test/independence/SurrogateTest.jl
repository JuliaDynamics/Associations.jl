
n = 100

# Pre-discretized data
likeit = rand(["yes", "no"], n)
food = rand(["veggies", "meat", "fish"], n)
service = rand(["netflix", "hbo"], n)
est = Contingency()
nsurr = 3

@test independence(SurrogateTest(MIShannon(), est; nsurr), food, likeit) isa SurrogateTestResult
@test independence(SurrogateTest(MIRenyiJizba(), est; nsurr), food, likeit) isa SurrogateTestResult
@test independence(SurrogateTest(MIRenyiSarbu(), est; nsurr), food, likeit) isa SurrogateTestResult
@test independence(SurrogateTest(MITsallisFuruichi(), est; nsurr), food, likeit) isa SurrogateTestResult
@test independence(SurrogateTest(MITsallisMartin(), est; nsurr), food, likeit) isa SurrogateTestResult

@test independence(SurrogateTest(CMIShannon(), est; nsurr), food, likeit, service) isa SurrogateTestResult
@test independence(SurrogateTest(CMIRenyiSarbu(), est; nsurr), food, likeit, service) isa SurrogateTestResult
@test independence(SurrogateTest(CMIRenyiJizba(), est; nsurr), food, likeit, service) isa SurrogateTestResult
