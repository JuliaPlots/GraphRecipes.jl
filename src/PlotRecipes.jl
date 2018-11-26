
module PlotRecipes

using LightGraphs
using NetworkLayout
using PlotUtils         # ColorGradient
using RecipesBase

using LinearAlgebra
using SparseArrays
using Statistics

include("utils.jl")
include("graph_layouts.jl")
include("graphs.jl")
include("misc.jl")

end # module
