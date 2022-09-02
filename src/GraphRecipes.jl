module GraphRecipes

using Graphs
using PlotUtils         # ColorGradient
using RecipesBase

using InteractiveUtils  # subtypes
using LinearAlgebra
using SparseArrays
using Statistics
using NaNMath
using GeometryTypes
using Interpolations

import NetworkLayout
import Graphs: getRNG

include("utils.jl")
include("graph_layouts.jl")
include("graphs.jl")
include("misc.jl")
include("trees.jl")

end
