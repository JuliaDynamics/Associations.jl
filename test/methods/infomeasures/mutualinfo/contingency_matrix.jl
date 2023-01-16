x = rand(["a", "b", "c"], 200)
y = rand(["hello", "yoyo", "heyhey"], 200)
c = contingency_matrix(x, y)

@test mutualinfo(MIShannon(), c) >= 0.0
@test mutualinfo(MITsallisFuruichi(), c) isa Real
@test mutualinfo(MITsallisMartin(), c) isa Real
@test mutualinfo(MIRenyiJizba(), c) isa Real
@test mutualinfo(MIRenyiSarbu(), c) isa Real
