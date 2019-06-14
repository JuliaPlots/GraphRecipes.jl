using Random
using VisualRegressionTests
using LinearAlgebra
using SparseArrays
using ImageMagick

istravis = "TRAVIS" ∈ keys(ENV)

default(show=false, reuse=true)

cd(@__DIR__)
cd("../assets")

@testset "README" begin
    n = 15
    Random.seed!(42)
    A = Float64[ rand() < 0.5 ? 0 : rand() for i=1:n, j=1:n]
    for i=1:n
        A[i, 1:i-1] = A[1:i-1, i]
    end
    x = rand(n)
    y = rand(n)
    z = rand(n)
    @plottest graphplot(A,
                        node_weights = 1:n,
                        marker = (:YlOrRd, :rect),
                        marker_z = 1:n,
                        markersize = 3,
                        names = 1:n,
                        linecolor = :darkgrey,
                        layout_kw = Dict(:x => x, :y => y),
                        ) "random_labelled_graph.png" popup=!istravis tol=0.03

    @plottest graphplot(A,
                        node_weights = 1:n,
                        markercolor = :darkgray,
                        dim = 3,
                        markersize = 5,
                        linecolor = :darkgrey,
                        linealpha = 0.5,
                        layout_kw = Dict(:x => x, :y => y, :z => z),
                        ) "random_3d_graph.png" popup=!istravis tol=0.02
    adjmat = Symmetric(sparse(rand(0:1,8,8)))

    adjmat = Symmetric(sparse(rand(0:1,8,8)))
    @plottest plot(
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
                  ) "arc_chord_diagrams.png" popup=!istravis tol=0.02

    code = :(
    function mysum(list)
        out = 0
        for value in list
            out += value
        end
        out
    end
    )

    @plottest plot(code, fontsize=5, shorten=0.2, axis_buffer=0.15, markersize=0) "AST_example.png" popup=!istravis tol=0.02
    @plottest plot(AbstractFloat, method=:tree, fontsize=4, markersize=0) "julia_type_tree.png" popup=!istravis tol=0.2
end
