using Random
using LinearAlgebra
using SparseArrays
using GraphRecipes
using Plots
using LightGraphs

cd(@__DIR__)
cd("../assets")

n = 15
Random.seed!(1)
A = Float64[ rand() < 0.5 ? 0 : rand() for i=1:n, j=1:n]
for i=1:n
    A[i, 1:i-1] = A[1:i-1, i]
end
x = rand(n)
y = rand(n)
z = rand(n)
graphplot(A,
          nodesize = 0.25,
          node_weights = 1:n,
          nodecolor = range(colorant"yellow", stop=colorant"red", length=n),
          names = 1:n,
          linecolor = :darkgrey,
          layout_kw = Dict(:x => x, :y => y),
          )
savefig("random_labelled_graph.png")

graphplot(A,
          node_weights = 1:n,
          markercolor = :darkgray,
          dim = 3,
          markersize = 5,
          linecolor = :darkgrey,
          linealpha = 0.5,
          layout_kw = Dict(:x => x, :y => y, :z => z),
          )
savefig("random_3d_graph.png")

g = wheel_graph(10)
graphplot(g, curves=false)
savefig("LightGraphs.png")

g = [0 1 1;
     0 0 1;
     0 1 0]

graphplot(g, names=1:3, curvature_scalar=0.1)
savefig("directed.png")

n = 8
g = wheel_digraph(n)
edgelabel_dict = Dict()
for i in 1:n
    for j in 1:n
        edgelabel_dict[(i, j)] = string("edge ", i, " to ", j)
    end
end

graphplot(g, names=1:n, edgelabel=edgelabel_dict, curves=false, nodeshape=:rect)

savefig("edgelabel.png")

g = [0 1 1;
     0 1 0;
     0 1 0]

graphplot(g)

savefig("selfedges.png")

graphplot([[1,1,2,2],[1,1,1],[1]], names="node_".*string.(1:3), nodeshape=:circle, self_edge_size=0.2)
savefig("multigraphs.png")

adjmat = Symmetric(sparse(rand(0:1,8,8)))
plot(
     graphplot(adjmat,
               method=:chorddiagram,
               names=[text(string(i), 8) for i in 1:8],
               linecolor=:black,
               fillcolor=:lightgray),

     graphplot(adjmat,
               method=:arcdiagram,
               markersize=3,
               linecolor=:black,
               markercolor=:black)
               )
savefig("arc_chord_diagrams.png")

code = :(
function mysum(list)
    out = 0
    for value in list
        out += value
    end
    out
end
)

default(size=(1000, 1000))
plot(code, fontsize=10, shorten=0.01, axis_buffer=0.15, nodeshape=:rect)
savefig("AST_example.png")

plot(AbstractFloat, method=:tree, fontsize=10, nodeshape=:ellipse)
savefig("julia_type_tree.png")


# `AbstractTrees` example
using AbstractTrees
AbstractTrees.children(d::Dict) = [p for p in d]
AbstractTrees.children(p::Pair) = AbstractTrees.children(p[2])
function AbstractTrees.printnode(io::IO, p::Pair)
    str = isempty(AbstractTrees.children(p[2])) ? string(p[1], ": ", p[2]) : string(p[1], ": ")
    print(io, str)
end

d = Dict(:a => 2,:d => Dict(:b => 4,:c => "Hello"),:e => 5.0)

plot(TreePlot(d), method=:tree, fontsize=10, nodeshape=:ellipse)
savefig("julia_dict_tree.png")
