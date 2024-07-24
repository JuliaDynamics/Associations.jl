
ests = [
    EntropyDecomposition(TEShannon(),  Lord()),
    EntropyDecomposition(CMIShannon(), Kraskov()),
]

@test occursin("EntropyDecomposition", repr(ests))