using Random
using VisualRegressionTests
using LinearAlgebra
using SparseArrays
using ImageMagick
using LightGraphs
using AbstractTrees

cd(@__DIR__)
include("functions.jl")

istravis = "TRAVIS" ∈ keys(ENV)

default(show=false, reuse=true)

cd("../assets")

@testset "FIGURES" begin
    @plottest random_labelled_graph() "random_labelled_graph.png" popup=!istravis tol=0.03

    @plottest random_3d_graph() "random_3d_graph.png" popup=!istravis tol=0.02

    @plottest light_graphs "light_graphs.png" popup=!istravis tol=0.02

    @plottest directed "directed.png" popup=!istravis tol=0.02

    @plottest marker_properties() "marker_properties.png" popup=!istravis tol=0.02

    @plottest edgelabel() "edgelabel.png" popup=!istravis tol=0.1

    @plottest selfedges() "selfedges.png" popup=!istravis tol=0.02

    @plottest multigraphs() "multigraphs.png" popup=!istravis tol=0.1

    @plottest arc_chord_diagrams() "arc_chord_diagrams.png" popup=!istravis tol=0.02

    @plottest ast_example() "ast_example.png" popup=!istravis tol=0.02

    @plottest julia_type_tree() "julia_type_tree.png" popup=!istravis tol=0.2

    @plottest julia_dict_tree() "julia_dict_tree.png" popup=!istravis tol=0.2

end
