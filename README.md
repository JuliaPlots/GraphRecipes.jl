# PlotRecipes

[![Build Status](https://travis-ci.org/JuliaPlots/PlotRecipes.jl.svg?branch=master)](https://travis-ci.org/JuliaPlots/PlotRecipes.jl)

### Primary author: Thomas Breloff (@tbreloff)

This repo maintains a collection of recipes for machine learning, graph analysis, finance, and more.  It uses the powerful machinery of [Plots](https://github.com/tbreloff/Plots.jl) and [RecipesBase](https://github.com/JuliaPlots/RecipesBase.jl) to turn simple transformations into flexible visualizations.

PlotRecipes also exports the recipes in [StatPlots.jl](https://github.com/JuliaPlots/StatPlots.jl), which is a collection of statistics recipes, including functionality for DataFrames and Distributions.

---

# Examples

---

## Graphs

#### Spectral

```julia
using PlotRecipes
n = 15
A = Float64[(rand()<0.5 ? 0 : rand()) for i=1:n,j=1:n]
for i=1:n
    A[i,1:i-1] = A[1:i-1,i]
end
node_weights = 1:n

graphplot(A,
    node_weights = 1:n,
    marker = (:heat, :rect),
    line = (3, 0.5, :blues),
    marker_z = 1:n,
    names = 1:n
)
```

![](https://cloud.githubusercontent.com/assets/933338/16093627/9da7b26a-330a-11e6-9733-9d28d5bab604.png)

```julia
graphplot(A,
    node_weights = 1:n,
    dim = 3,
    line = (3, 0.5, :blues),
    marker_z = 1:n
)
```

![](https://cloud.githubusercontent.com/assets/933338/16094180/0dd2edf0-330d-11e6-8596-d12b0b8d5393.png)

#### Arc and chord diagrams

```julia
using PlotRecipes

adjmat = Symmetric(sparse(rand(0:1,8,8)));

plot(
    graphplot(adjmat, method=:chorddiagram),
    graphplot(adjmat, method=:arcdiagram, markersize=3)
    )
```  
![arc and chord diagrams](https://user-images.githubusercontent.com/2822757/27743452-5511e5e2-5dbc-11e7-895e-dfa753a84efc.png)  

#### Fun with algos. Visualizing a stress-driven layout algorithm

![](https://cloud.githubusercontent.com/assets/933338/16698919/ee1f9e76-451e-11e6-8936-881551f120dd.gif)

---

#### Julia code -- AST

```julia
using PlotRecipes
pyplot(ma=0.8,lc=:white,mc=:white,size=(1000,800))
theme(:dark)

code = :(
function transform!{T}(act::Activation{:softmax,T})
    val = forward!(act.input)
    out = act.output.val
    for i=1:act.n
        out[i] = exp(val[i])
    end
    s = one(T) / sum(out)
    for i=1:act.n
        out[i] *= s
    end
    out
end
)

plot(code, fontsize=11, shorten=0.2, axis_buffer=0.05)
```

![](https://cloud.githubusercontent.com/assets/933338/20402948/cb618014-accc-11e6-969a-28e738a8bea0.png)

---

#### Julia Type Trees

```julia
using PlotRecipes, Learn
pyplot(size=(800,500))
theme(:dark)
plot(Learnable, method=:tree)
```

![](https://cloud.githubusercontent.com/assets/933338/20758853/2420f72c-b6e9-11e6-82dd-4e62a679b3cb.png)

---

## Maps and Shapefile.jl

```julia
using PlotRecipes, Shapefile
dir = "https://github.com/nvkelso/natural-earth-vector/raw/master/110m_physical/"
fn = "ne_110m_land.shp"
run(`wget $dir/$fn -P /tmp/`)
shp = open("/tmp/$fn") do fd
    read(fd, Shapefile.Handle)
end
shapeplot(shp.shapes, c=:grey)
```

![](https://cloud.githubusercontent.com/assets/933338/16770876/83dea362-481c-11e6-9943-bb77148be5b8.png)

---

#### Portfolio Composition Maps

```julia
using PlotRecipes
tickers = ["IBM", "Google", "Apple", "Intel"]
N = 10
D = length(tickers)
weights = rand(N,D)
weights ./= sum(weights, 2)
returns = sort!((1:N) + D*randn(N))

portfoliocomposition(weights, returns, labels = tickers')
```

![](https://github.com/JuliaPlots/PlotReferenceImages.jl/blob/master/PlotRecipes/portfoliocomposition.png)

---

## StatPlots.jl


---

#### corrplot

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
using PlotRecipes, Distributions
n = 1000
x = rand(Gamma(2), n)
y = -0.5x + randn(n)

marginalhist(x, y, fc=:plasma, bins=40)
```

![](https://github.com/JuliaPlots/PlotReferenceImages.jl/blob/master/PlotRecipes/marginalhist.png)
