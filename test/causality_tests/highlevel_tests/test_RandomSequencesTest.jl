import UncertainData: RandomSequences 

rtest = RandomSequencesTest(CrossMappingTest(), RandomSequences(10, 10))
@test rtest isa RandomSequencesTest