
__precompile__()

module PlotRecipes

using Reexport
@reexport using Plots
@reexport using StatPlots

include("graphs.jl")
# include("hists.jl")  # TODO: remove when StatPlots is ready
# include("stats.jl")  # TODO: remove when StatPlots is ready
include("finance.jl")
include("misc.jl")

end # module
