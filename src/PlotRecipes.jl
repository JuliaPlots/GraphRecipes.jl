
module PlotRecipes

using Reexport
@reexport using Plots
import Plots: Plot, isnothing
import NetworkLayout

include("utils.jl")
include("graphs.jl")
include("finance.jl")
include("misc.jl")

@deprecate arcdiagram graphplot
@deprecate chorddiagram graphplot

end # module
