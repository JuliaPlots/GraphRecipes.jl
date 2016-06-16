
# see: http://stackoverflow.com/a/37732384/5075246

@userplot PortfolioComposition

# this shows the shifting composition of a basket of something over a variable
# - "returns" are the dependent variable
# - "weights" are a matrix where the ith column is the composition for returns[i]
# - since each polygon is its own series, you can assign labels easily
@recipe function f(pc::PortfolioComposition)
    weights, returns = pc.args
    n = length(returns)
    weights = cumsum(weights,2)
    seriestype := :shape

 #    # Set the axis ticks for non-numeric data
 #    try
 #    	if !(eltype(returns) <: Number)
 #    		ticksym = Plots.isvertical(d) ? :yticks : :xticks
 #    		d[ticksym] = ((1:n)-0.5, returns)
	# 	end
	# end

	# create a filled polygon for each item
    for c=1:size(weights,2)
        sx = vcat(weights[:,c], c==1 ? zeros(n) : reverse(weights[:,c-1]))
        sy = vcat(returns, reverse(returns))
        @series Plots.isvertical(d) ? (sx, sy) : (sy, sx)
        # @series Plots.isvertical(d) ? Shape(sx, sy) : Shape(sy, sx)
    end
end

