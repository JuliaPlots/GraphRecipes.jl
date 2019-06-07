using GraphRecipes
using Plots
using Test

@testset "utils.jl" begin

    @test GraphRecipes.directed_curve(0., 1., 0., 1.) == GraphRecipes.directed_curve(0, 1, 0, 1)

    @testset "Functions from Plots.jl" begin

        @test GraphRecipes.isnothing(nothing) == Plots.isnothing(nothing)
        @test GraphRecipes.isnothing(missing) == Plots.isnothing(missing)
        @test GraphRecipes.isnothing(NaN) == Plots.isnothing(NaN)
        @test GraphRecipes.isnothing(0) == Plots.isnothing(0)
        @test GraphRecipes.isnothing(1) == Plots.isnothing(1)
        @test GraphRecipes.isnothing(0.0) == Plots.isnothing(0.0)
        @test GraphRecipes.isnothing(1.0) == Plots.isnothing(1.0)

        for (s, e) in [ (rand(), rand()) for i in 1:100 ]
            @test GraphRecipes.partialcircle(s, e) == Plots.partialcircle(s, e)
        end

    end

end

@testset "spectral" begin
    adjmat = [0  1  0  1;
              0  0  1  1;
              0  0  0  0;
              1  0  1  0]

    x, y, z = GraphRecipes.spectral_graph(adjmat)
    @test x == [-0.9999999914991978;
                -0.5000000055998736;
                 0.35146814988501696;
                 0.9999999933743652]
    @test y == [-0.9999999914991975;
                 0.9226039113558594;
                 0.9999999524609413;
                -0.5999999960246193]
    @test z == [-0.9999999914991974;
                -0.8363019304325268;
                -0.23639887033225468;
                -0.600000023554825]
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
