using PlotRecipes
using Base.Test

# write your own tests here
@test 1 == 1

# TODO: actual tests

using Distributions
n = 1000
x = rand(Gamma(2), n)
y = -0.5x + randn(n)
marginalhist(x, y)