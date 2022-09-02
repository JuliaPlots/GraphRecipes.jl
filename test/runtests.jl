using VisualRegressionTests
using AbstractTrees
using LinearAlgebra
using GraphRecipes
using SparseArrays
using ImageMagick
using Graphs
using Plots
using Test
using Gtk  # for popup

import Graphs: getRNG

isci() = get(ENV, "CI", "false") == "true"
itol(tol = nothing) = something(tol, isci() ? 1e-3 : 1e-5)

include("functions.jl")
include("parse_readme.jl")

default(show=false, reuse=true)

cd(joinpath(@__DIR__, "..", "assets")) do
    @testset "FIGURES" begin
        @plottest random_labelled_graph() "random_labelled_graph.png" popup=!isci() tol=itol()

        @plottest random_3d_graph() "random_3d_graph.png" popup=!isci() tol=itol()

        @plottest light_graphs() "light_graphs.png" popup=!isci() tol=itol()

        @plottest directed() "directed.png" popup=!isci() tol=itol()

        @plottest marker_properties() "marker_properties.png" popup=!isci() tol=itol()

        @plottest edgelabel() "edgelabel.png" popup=!isci() tol=itol()

        @plottest selfedges() "selfedges.png" popup=!isci() tol=itol()

        @plottest multigraphs() "multigraphs.png" popup=!isci() tol=itol()

        @plottest arc_chord_diagrams() "arc_chord_diagrams.png" popup=!isci() tol=itol()

        @plottest ast_example() "ast_example.png" popup=!isci() tol=itol()

        @plottest julia_type_tree() "julia_type_tree.png" popup=!isci() tol=itol()

        @plottest julia_dict_tree() "julia_dict_tree.png" popup=!isci() tol=itol()

        @plottest funky_edge_and_marker_args() "funky_edge_and_marker_args.png" popup=!isci() tol=itol()

        @plottest custom_nodeshapes_single() "custom_nodeshapes_single.png" popup=!isci() tol=itol()

        @plottest custom_nodeshapes_various() "custom_nodeshapes_various.png" popup=!isci() tol=itol()
    end

    @testset "README" begin
        @plottest julia_logo_pun() "readme_julia_logo_pun.png" popup=!isci tol=itol()
    end
end

@testset "utils.jl" begin
    @test GraphRecipes.directed_curve(0., 1., 0., 1.) == GraphRecipes.directed_curve(0, 1, 0, 1)

    @testset "Functions from Plots.jl" begin

        @test GraphRecipes.isnothing(nothing) == Plots.isnothing(nothing)
        @test GraphRecipes.isnothing(missing) == Plots.isnothing(missing)
        @test GraphRecipes.isnothing(NaN) == Plots.isnothing(NaN)
        @test GraphRecipes.isnothing(0) == Plots.isnothing(0)
        @test GraphRecipes.isnothing(1) == Plots.isnothing(1)
        @test GraphRecipes.isnothing(0.0) == Plots.isnothing(0.0)
        @test GraphRecipes.isnothing(1.0) == Plots.isnothing(1.0)

        rng = getRNG(123)
        for (s, e) in [ (rand(rng), rand(rng)) for i in 1:100 ]
            @test GraphRecipes.partialcircle(s, e) == Plots.partialcircle(s, e)
        end

        @testset "nearest_intersection" begin
            @test GraphRecipes.nearest_intersection(0, 0, 3, 3, [(1,0), (0,1)]) == (0,0,0.5,0.5)
            @test GraphRecipes.nearest_intersection(1, 2, 1, 2, []) == (1, 2, 1, 2)
        end

        @testset "unoccupied_angle" begin
            @test GraphRecipes.unoccupied_angle(1, 1, [1, 1, 1, 1], [2, 0, 3, -1]) == 2pi
        end

        @testset "islabel" begin
            @test GraphRecipes.islabel("hi")
            @test GraphRecipes.islabel(1)
            @test !GraphRecipes.islabel(missing)
            @test !GraphRecipes.islabel(NaN)
            @test !GraphRecipes.islabel(false)
            @test !GraphRecipes.islabel("")
        end

        @testset "control_point" begin
            @test GraphRecipes.control_point(0, 0, 6, 0, 4) == (4,3)
        end

        # TODO: Actually test that the aliases produce the same plots, rather than just
        # checking that they don't error. Also, test all of the different aliases.
        @testset "Aliases" begin
            A = [1 0 1 0;0 0 1 1;1 1 1 1;0 0 1 1]
            graphplot(A, markercolor=:red, markershape=:rect, markersize=0.5)
            graphplot(A, nodeweights=1:4)
            graphplot(A, curvaturescalar=0)
            graphplot(A, el=Dict((1,2)=>""), elb=true)
            graphplot(A, ew=(s,d,w)->3)
            graphplot(A, ses=0.5)
        end
    end
end

# -----------------------------------------
# marginalhist

# rng = getRNG(123)
# using Distributions
# n = 1000
# x = rand(rng, Gamma(2), n)
# y = -0.5x + randn(rng, n)
# marginalhist(x, y)

# -----------------------------------------
# portfolio composition map

# # fake data
# rng = getRNG(123)
# tickers = ["IBM", "Google", "Apple", "Intel"]
# N = 10
# D = length(tickers)
# weights = rand(N,D)
# weights ./= sum(weights, 2)
# returns = sort!((1:N) + D*randn(rng, N))

# # plot it
# portfoliocomposition(weights, returns, labels = tickers')

# -----------------------------------------
#
