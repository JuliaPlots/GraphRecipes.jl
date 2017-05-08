
module PlotRecipes

using Reexport
@reexport using StatPlots
@reexport using Plots
import NetworkLayout

include("utils.jl")
include("graphs.jl")
include("finance.jl")
include("misc.jl")

end # module
