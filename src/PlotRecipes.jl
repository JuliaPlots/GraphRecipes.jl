
module PlotRecipes

using Reexport
@reexport using Plots
import Plots: Plot
@reexport using StatPlots
import NetworkLayout

include("graphs.jl")
include("finance.jl")
include("misc.jl")

end # module
