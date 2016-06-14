# PlotRecipes

### Primary author: Thomas Breloff (@tbreloff)

This repo maintains a collection of recipes for statistics, machine learning, graph analysis, finance, and more.  It uses the powerful machinery of [Plots](https://github.com/tbreloff/Plots.jl) and [RecipesBase](https://github.com/JuliaPlots/RecipesBase.jl) to turn simple transformations into flexible visualizations.

---

# Examples

---

#### Stats

```julia
using PlotRecipes
M = randn(1000,4)
M[:,2] += 0.8sqrt(abs(M[:,1])) - 0.5M[:,3] + 5
M[:,3] -= 0.7M[:,1].^2 + 2

corrplot(M, label = ["x$i" for i=1:4])
```

![](https://cloud.githubusercontent.com/assets/933338/16030833/3c84e6bc-31c3-11e6-9a04-4cee531440a4.png)

--

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
