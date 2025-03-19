using Graphs
using GraphPlot, Compose
g = SimpleGraph(2)
add_edge!(g, 1, 2)
gplot(g,locs_x=[1,2],locs_y=[1,2])