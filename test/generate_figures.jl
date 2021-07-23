using StableRNGs
using LinearAlgebra
using SparseArrays
using GraphRecipes
using Plots
using LightGraphs
using AbstractTrees

cd(@__DIR__)
include("parse_readme.jl")
include("functions.jl")
cd("../assets")

julia_logo_pun()
savefig("readme_julia_logo_pun.png")

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

funky_edge_and_marker_args()
savefig("funky_edge_and_marker_args.png")

custom_nodeshapes_single()
savefig("custom_nodeshapes_single.png")

custom_nodeshapes_various()
savefig("custom_nodeshapes_various.png")
