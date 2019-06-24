using Random
using LinearAlgebra
using SparseArrays
using GraphRecipes
using Plots

cd(@__DIR__)
cd("../assets")

n = 15
Random.seed!(42)
A = Float64[ rand() < 0.5 ? 0 : rand() for i=1:n, j=1:n]
for i=1:n
    A[i, 1:i-1] = A[1:i-1, i]
end
x = rand(n)
y = rand(n)
z = rand(n)
graphplot(A,
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

adjmat = Symmetric(sparse(rand(0:1,8,8)))

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

plot(code, fontsize=5, shorten=0.01, axis_buffer=0.15, nodeshape=:rect)
savefig("AST_example.png")

plot(AbstractFloat, method=:tree, fontsize=4, nodeshape=:rectangle)
savefig("julia_type_tree.png")
