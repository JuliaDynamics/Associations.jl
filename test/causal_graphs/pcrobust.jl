using Graphs
x, y, z = randn(100), randn(100), randn(100)
g = infer_graph(PCRobust())
