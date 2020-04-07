using Random
using LinearAlgebra
using SparseArrays
using GraphRecipes
using Plots
using LightGraphs
using AbstractTrees

cd(@__DIR__)
include("functions.jl")
cd("../assets")

random_labelled_graph()[1]
savefig("random_labelled_graph.png")

random_3d_graph()
savefig("random_3d_graph.png")

light_graphs()
savefig("light_graphs.png")

directed()
savefig("directed.png")

edgelabel()
savefig("edgelabel.png")

selfedges()
savefig("selfedges.png")

multigraphs()
savefig("multigraphs.png")

arc_chord_diagrams()
savefig("arc_chord_diagrams.png")

marker_properties()
savefig("marker_properties.png")

ast_example()
savefig("ast_example.png")

julia_type_tree()
savefig("julia_type_tree.png")

julia_dict_tree()
savefig("julia_dict_tree.png")
