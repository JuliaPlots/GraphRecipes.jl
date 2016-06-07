# PlotRecipes

See [Plots](https://github.com/tbreloff/Plots.jl) and [RecipesBase](https://github.com/JuliaPlots/RecipesBase.jl) for more info.  Examples:

```julia
using Distributions, PlotRecipes
pyplot()
n = 1000
x = rand(Gamma(2), n)
y = -0.5x + randn(n)
marginalhist(x, y, fc=:plasma, bins=40)
```

![](https://github.com/JuliaPlots/PlotReferenceImages.jl/blob/master/PlotRecipes/marginalhist.png)
