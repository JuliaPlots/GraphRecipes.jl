using PlotRecipes
using Test

@testset "utils.jl" begin

    @test directed_curve(0.,1.,0.,1.) == directed_curve(0,1,0,1)

    
end

# -----------------------------------------
# marginalhist

# using Distributions
# n = 1000
# x = rand(Gamma(2), n)
# y = -0.5x + randn(n)
# marginalhist(x, y)

# -----------------------------------------
# portfolio composition map

# # fake data
# tickers = ["IBM", "Google", "Apple", "Intel"]
# N = 10
# D = length(tickers)
# weights = rand(N,D)
# weights ./= sum(weights, 2)
# returns = sort!((1:N) + D*randn(N))

# # plot it
# portfoliocomposition(weights, returns, labels = tickers')

# -----------------------------------------
#
