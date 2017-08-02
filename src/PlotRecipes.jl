
module PlotRecipes

using Reexport
@reexport using Plots
import Plots: Plot, isnothing
@reexport using StatPlots
import NetworkLayout

include("utils.jl")
include("graphs.jl")
include("finance.jl")
include("misc.jl")

end # module
