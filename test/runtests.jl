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
    @testset "edge_label_exsists" begin
        @test GraphRecipes.edge_label_exsists(Dict((1,2) => "hello"), (1,3)) == false
        @test GraphRecipes.edge_label_exsists(Dict((1,3) => "hi"), (1,3)) == true
        @test GraphRecipes.edge_label_exsists(["hi",2,4,missing,false,nothing,1,2,3], (2,1)) == true
        @test GraphRecipes.edge_label_exsists(["hi",2,4,missing,false,nothing,1,2,3], (2,2)) == false
        @test GraphRecipes.edge_label_exsists(["hi",2,4,missing,false,nothing,1,2,3], (1,2)) == false
        @test GraphRecipes.edge_label_exsists(["hi",2,4,missing,false,nothing,1,2,3], (3,2)) == false
        @test GraphRecipes.edge_label_exsists(["hi",2,4,missing,false,NaN,1,2,3], (3,2)) == false
        @test GraphRecipes.edge_label_exsists(["hi",2,4,missing,false,"",1,2,3], (3,2)) == false
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
