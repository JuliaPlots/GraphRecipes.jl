# GraphRecipes
The repository formerly known as PlotRecipes

[![Build Status](https://travis-ci.org/JuliaPlots/GraphRecipes.jl.svg?branch=master)](https://travis-ci.org/JuliaPlots/GraphRecipes.jl)
[Documentation](http://docs.juliaplots.org/latest/graphrecipes/introduction/)

## Summary
In this repository, a graph is a network of connected nodes (although sometimes people use the same word to refer to a plot). If you want to do plotting, then use [Plots.jl](https://github.com/JuliaPlots/Plots.jl).

For a given graph, there are many legitimate ways to display and visualize that graph. However, some graph layouts will convey the structure of the underlying graph much more clearly than other layouts. GraphRecipes provides many options for producing graph layouts including  (un)directed graphs, tree graphs and arc/chord diagrams. For each layout type the `graphplot` function will try to create a default layout that optimizes visual clarity. However, the user can tweak the default layout through a large number of powerful keyword arguments, see the [Documentation](http://docs.juliaplots.org/latest/graphrecipes/introduction/) for more details and some examples.

## Installation
```julia
]add GraphRecipes
```

## An example
```julia
using GraphRecipes
using Plots

g = [0 1 1;
     1 0 1;
     1 1 0]

graphplot(g,
          x=[0,-1/tan(π/3),1/tan(π/3)], y=[1,0,0],
          nodeshape=:circle, nodesize=1.1,
          axis_buffer=0.6,
          curves=false,
          color=:black,
          nodecolor=[colorant"#389826",colorant"#CB3C33",colorant"#9558B2"],
          linewidth=10)
```
![juliagraph](https://user-images.githubusercontent.com/8610352/76132303-27016900-6077-11ea-8fb1-d613c4ddda59.png)

Original author: Thomas Breloff (@tbreloff)

This repo maintains a collection of recipes for graph analysis, and is a reduced and refactored version of the previous PlotRecipes. It uses the powerful machinery of [Plots](https://github.com/JuliPlots/Plots.jl) and [RecipesBase](https://github.com/JuliaPlots/RecipesBase.jl) to turn simple transformations into flexible visualizations.
