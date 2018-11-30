using PlotRecipes
using Plots
using Test

@testset "utils.jl" begin

    @test PlotRecipes.directed_curve(0., 1., 0., 1.) == PlotRecipes.directed_curve(0, 1, 0, 1)

    @testset "Functions from Plots.jl" begin

        @test PlotRecipes.isnothing(nothing) == Plots.isnothing(nothing)
        @test PlotRecipes.isnothing(missing) == Plots.isnothing(missing)
        @test PlotRecipes.isnothing(NaN) == Plots.isnothing(NaN)
        @test PlotRecipes.isnothing(0) == Plots.isnothing(0)
        @test PlotRecipes.isnothing(1) == Plots.isnothing(1)
        @test PlotRecipes.isnothing(0.0) == Plots.isnothing(0.0)
        @test PlotRecipes.isnothing(1.0) == Plots.isnothing(1.0)

        for (s, e) in [ (rand(), rand()) for i in 1:100 ]
            @test PlotRecipes.partialcircle(s, e) == Plots.partialcircle(s, e)
        end

    end

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
