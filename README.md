# PlotRecipes

### Primary author: Thomas Breloff (@tbreloff)

This repo maintains a collection of recipes for statistics, machine learning, graph analysis, finance, and more.  It uses the powerful machinery of [Plots](https://github.com/tbreloff/Plots.jl) and [RecipesBase](https://github.com/JuliaPlots/RecipesBase.jl) to turn simple transformations into flexible visualizations.

---

# Examples

---

#### Graphs



---

#### Marginal Histograms

```julia
using Distributions, PlotRecipes
pyplot()
n = 1000
x = rand(Gamma(2), n)
y = -0.5x + randn(n)
marginalhist(x, y, fc=:plasma, bins=40)
```

![](https://github.com/JuliaPlots/PlotReferenceImages.jl/blob/master/PlotRecipes/marginalhist.png)

---

#### Portfolio Composition Maps

```julia
using PlotRecipes

# fake data
tickers = ["IBM", "Google", "Apple", "Intel"]
N = 10
D = length(tickers)
weights = rand(N,D)
weights ./= sum(weights, 2)
returns = sort!((1:N) + D*randn(N))

# plot it
portfoliocomposition(weights, returns, labels = tickers')
```

![](https://github.com/JuliaPlots/PlotReferenceImages.jl/blob/master/PlotRecipes/portfoliocomposition.png)

---
