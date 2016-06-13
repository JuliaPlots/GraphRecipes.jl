module PlotRecipes

using Reexport
@reexport using Plots

include("graphs.jl")
include("hists.jl")
include("stats.jl")

end # module
