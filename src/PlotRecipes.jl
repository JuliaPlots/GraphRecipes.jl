
__precompile__()

module PlotRecipes

using Reexport
@reexport using Plots
@reexport using StatPlots

include("graphs.jl")
include("finance.jl")
include("misc.jl")

end # module
