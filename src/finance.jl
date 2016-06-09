
# see: http://stackoverflow.com/a/37732384/5075246

@userplot PortfolioComposition

@recipe function f(pc::PortfolioComposition)
    weights, returns = pc.args
    weights = cumsum(weights,2)
    seriestype := :shape
    for c=1:size(weights,2)
        sx = vcat(weights[:,c], c==1 ? zeros(length(returns)) : reverse(weights[:,c-1]))
        sy = vcat(returns, reverse(returns))
        @series Shape(sx, sy)
    end
end

