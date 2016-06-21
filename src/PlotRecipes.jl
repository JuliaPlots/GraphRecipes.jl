module PlotRecipes

using Reexport
@reexport using Plots

include("graphs.jl")
include("hists.jl")
include("stats.jl")
include("finance.jl")
include("misc.jl")

end # module
