# @userplot MarginalHist
@shorthands marginalhist

# @recipe function f(h::MarginalHist)
#     if length(h.args) != 2 || !(typeof(h.args[1]) <: AbstractVector) || !(typeof(h.args[2]) <: AbstractVector)
#         error("Marginal Histograms should be given two vectors.  Got: $(typeof(h.args))")
#     end
#     x, y = h.args

@recipe function f(::Type{Val{:marginalhist}}, plt::Plot; density = false)
    x, y = d[:x], d[:y]
    delete!(d, :density)
    
    # set up the subplots
    legend --> false
    link := :both
    grid --> false
    layout --> @layout [
        tophist           _
        hist2d{0.9w,0.9h} righthist
    ]
    
    # main histogram2d
    @series begin
        seriestype := :histogram2d
        subplot := 2
    end
    
    # these are common to both marginal histograms... borrow the first color from the fill gradient
    ticks --> nothing
    foreground_color_subplot := RGBA(0,0,0,0)
    fillcolor := getColor(colorscheme(get(d, :fillcolor, Plots.default_gradient())))
    fillalpha --> 0.3
    linealpha --> 0.3

    if density
        trim := true
        seriestype := :density
    else
        seriestype := :histogram
    end
    
    # upper histogram
    @series begin
        subplot := 1
        y := x
    end
    
    # right histogram
    @series begin
        orientation := :h
        subplot := 3
        y := y
    end
end

# # now you can plot like:
# marginalhist(rand(1000), rand(1000))
