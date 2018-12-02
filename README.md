# GraphRecipes
The repository formerly know as PlotRecipes

[![Build Status](https://travis-ci.org/JuliaPlots/GraphRecipes.jl.svg?branch=master)](https://travis-ci.org/JuliaPlots/GraphRecipes.jl)

### Primary author: Thomas Breloff (@tbreloff)

This repo maintains a collection of recipes for graph analysis, and is a reduced and refactored version of the previous PlotRecipes. It uses the powerful machinery of [Plots](https://github.com/tbreloff/Plots.jl) and [RecipesBase](https://github.com/JuliaPlots/RecipesBase.jl) to turn simple transformations into flexible visualizations.

# Examples

#### Spectral

```julia
using GraphRecipes
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

![graph_one](https://user-images.githubusercontent.com/2822757/49309894-072adf00-f4dd-11e8-8e4f-0d6c4d3de77c.png)

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

![graph_two](https://user-images.githubusercontent.com/2822757/49309891-02fec180-f4dd-11e8-999a-9a4d68e9e0a9.png)

#### Arc and chord diagrams

```julia
using LinearAlgebra
using SparseArrays
using GraphRecipes
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
![graph_three](https://user-images.githubusercontent.com/2822757/49309879-f9755980-f4dc-11e8-99c6-545f0e44f118.png)

#### Julia code -- AST

```julia
using GraphRecipes
using Plots
pyplot(ma=0.8,lc=:white,mc=:white,size=(800,600))
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

![graph_four](https://user-images.githubusercontent.com/2822757/49310588-faa78600-f4de-11e8-95cf-4587d0ba1077.png)

#### Julia Type Trees

```julia
using GraphRecipes
using Plots

pyplot(size=(800,600))
theme(:dark)

plot(Integer, method=:tree, fontsize=4)

```
![graph_five](https://user-images.githubusercontent.com/2822757/49309857-e3679900-f4dc-11e8-8b57-f878a6d9cb5e.png)
