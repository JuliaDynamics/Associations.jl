using DelayEmbeddings
#x, y, z, w = 1:30, 10:10:300, 100:100:3000, 1000:1000:30000
#data = Dataset(x, y, z, w)

# This actually works!
# x, y, z = 1:10, 1:10, 1:10
# js = (1, 2, 3)
# τs = (0, -3, 2)
# genembed(Dataset(x, y, z), τs, js)

x, y, z, w = rand(1000), rand(1000), rand(1000), rand(1000)
# This works now
targets, sources = preprocess_pa([x, y], [z, w]; τmax = 2, ηmax = 1)
est = VisitationFrequency(RectangularBinning(3))
te_aut(targets, sources, est)
