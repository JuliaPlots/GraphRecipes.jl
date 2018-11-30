# PlotRecipes

[![Build Status](https://travis-ci.org/JuliaPlots/PlotRecipes.jl.svg?branch=master)](https://travis-ci.org/JuliaPlots/PlotRecipes.jl)

### Primary author: Thomas Breloff (@tbreloff)

This repo maintains a collection of recipes for graph analysis. It uses the powerful machinery of [Plots](https://github.com/tbreloff/Plots.jl) and [RecipesBase](https://github.com/JuliaPlots/RecipesBase.jl) to turn simple transformations into flexible visualizations.

---

# Examples

---

## Graphs

#### Spectral

```julia
using PlotRecipes
using Plots

const n = 15
const A = Float64[ rand() < 0.5 ? 0 : rand() for i=1:n, j=1:n]
for i=1:n
    A[i, 1:i-1] = A[1:i-1, i]
end

graphplot(A,
          node_weights = 1:n,
          marker = (:YlOrRd, :rect),
          marker_z = 1:n,
          markersize = 3,
          names = 1:n,
          linecolor = :darkgrey,
       )

```

![](https://cloud.githubusercontent.com/assets/933338/16093627/9da7b26a-330a-11e6-9733-9d28d5bab604.png)

```julia
graphplot(A,
           node_weights = 1:n,
           markercolor = :darkgray,
           dim = 3,
           markersize = 5,
           linecolor = :darkgrey,
           linealpha = 0.5
       )

```

![](https://cloud.githubusercontent.com/assets/933338/16094180/0dd2edf0-330d-11e6-8596-d12b0b8d5393.png)

#### Arc and chord diagrams

```julia
using LinearAlgebra
using SparseArrays
using PlotRecipes
using Plots

adjmat = Symmetric(sparse(rand(0:1,8,8)));

plot(
    graphplot(adjmat,
              method=:chorddiagram,
              names=[text(string(i), 8) for i in 1:8],
              linecolor=:black,
              fillcolor=:lightgray),

    graphplot(adjmat,
              method=:arcdiagram,
              markersize=3,
              linecolor=:black,
              markercolor=:black)
    )

```
![arc and chord diagrams](https://user-images.githubusercontent.com/2822757/27743452-5511e5e2-5dbc-11e7-895e-dfa753a84efc.png)

#### Fun with algos. Visualizing a stress-driven layout algorithm

![](https://cloud.githubusercontent.com/assets/933338/16698919/ee1f9e76-451e-11e6-8936-881551f120dd.gif)

---

#### Julia code -- AST

```julia
using PlotRecipes
using Plots
pyplot(ma=0.8,lc=:white,mc=:white,size=(1000,800))
theme(:dark)

code = :(
function mysum(list)
    out = 0
    for value in list
        out += value
    end
    out
end
)

plot(code, fontsize=5, shorten=0.2, axis_buffer=0.05)

```

![](https://cloud.githubusercontent.com/assets/933338/20402948/cb618014-accc-11e6-969a-28e738a8bea0.png)

---

#### Julia Type Trees

```julia
using PlotRecipes
using Plots

pyplot(size=(800,500))
theme(:dark)

plot(Number, method=:tree)

```

![](https://cloud.githubusercontent.com/assets/933338/20758853/2420f72c-b6e9-11e6-82dd-4e62a679b3cb.png)
